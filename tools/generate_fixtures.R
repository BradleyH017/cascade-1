#!/usr/bin/env Rscript
#
# Generate validation fixtures for CASCADE testing
#
# Samples genes, peaks, and variants from production raw data using
# stratified sampling based on raw-data predicates (not package output).
# Extracts minimal inputs and runs the oracle to derive expected outputs.
#
# Usage:
#   Rscript tools/generate_fixtures.R <config_json> [output_dir]
#
# Output: .rds fixture files in tests/testthat/fixtures/

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# Null-coalescing operator (base R doesn't have %||% before R 4.4)
`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript generate_fixtures.R <config_json> [output_dir]")
}

config_file <- args[1]
output_dir <- if (length(args) >= 2) args[2] else "tests/testthat/fixtures"

config <- fromJSON(config_file, simplifyVector = FALSE)
cat("Generating fixtures from config:", config_file, "\n")
cat("Output directory:", output_dir, "\n\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Source the oracle (no cascade imports)
oracle_path <- file.path(dirname(dirname(output_dir)), "tests/testthat/helper-oracle.R")
if (!file.exists(oracle_path)) {
  # Try relative to package root
  oracle_path <- "tests/testthat/helper-oracle.R"
}
source(oracle_path)

set.seed(42)

# ============================================================================
# Helper: construct file path from pattern
# ============================================================================
construct_path <- function(pattern, cell_type = NULL, chr = NULL) {
  path <- pattern
  if (!is.null(cell_type)) path <- gsub("\\{CELL_TYPE\\}", cell_type, path)
  if (!is.null(chr)) path <- gsub("\\{CHR\\}", chr, path)
  path
}

# ============================================================================
# Config extraction
# ============================================================================
cell_types <- unlist(config$cell_types)
l1_cts <- cell_types[grepl("l1\\.", cell_types)]
l1_cts_no_pbmc <- l1_cts[!grepl("PBMC", l1_cts)]
chromosomes <- unlist(config$chromosomes)
eqtl_acat_pattern <- config$file_patterns$eqtl_acat
caqtl_acat_pattern <- config$file_patterns$caqtl_acat
eqtl_lfsr_file <- config$file_patterns$eqtl_lfsr
caqtl_lfsr_file <- config$file_patterns$caqtl_lfsr

# Column mappings
eqtl_acat_id_col <- config$column_mapping$eqtl_acat$feature_id %||% "phenotype_id"
eqtl_acat_q_col <- config$column_mapping$eqtl_acat$q_value %||% "ACAT_q"
caqtl_acat_id_col <- config$column_mapping$caqtl_acat$feature_id %||% "phenotype_id"
caqtl_acat_q_col <- config$column_mapping$caqtl_acat$q_value %||% "qval"

sig_threshold <- config$parameters$lfsr_sig_threshold %||% 0.05

cat("Cell types:", length(cell_types), "(", length(l1_cts), "L1 )\n")
cat("Chromosomes:", length(chromosomes), "\n")
cat("Sig threshold:", sig_threshold, "\n\n")

# ============================================================================
# GENE FIXTURES (N=20)
# ============================================================================
cat("=== Generating gene fixtures ===\n")

# Load eQTL ACAT q-values per L1 cell type
gene_acat_list <- list()
for (ct in l1_cts) {
  f <- construct_path(eqtl_acat_pattern, cell_type = ct)
  if (file.exists(f)) {
    dt <- fread(f, select = c(eqtl_acat_id_col, eqtl_acat_q_col))
    setnames(dt, c(eqtl_acat_id_col, eqtl_acat_q_col), c("gene_id", "q_value"))
    dt[, cell_type := ct]
    gene_acat_list[[ct]] <- dt
  }
}
gene_acat <- rbindlist(gene_acat_list)

# Compute per-gene significance profile
gene_sig <- gene_acat[q_value < sig_threshold, .(sig_cts = list(cell_type)), by = gene_id]
gene_all <- gene_acat[, .(tested_cts = list(unique(cell_type)), n_tested = uniqueN(cell_type)), by = gene_id]
gene_prof <- merge(gene_all, gene_sig, by = "gene_id", all.x = TRUE)
gene_prof[is.na(sig_cts), sig_cts := list(list(character(0)))]

# Map to L1 and determine lineage presence
myeloid_set <- ORACLE_MYELOID
lymphoid_set <- ORACLE_LYMPHOID

gene_prof[, sig_l1 := lapply(sig_cts, function(x) unique(oracle_map_to_l1(unlist(x))))]
gene_prof[, n_sig_l1 := sapply(sig_l1, function(x) length(setdiff(x, ORACLE_PBMC)))]
gene_prof[, has_myeloid := sapply(sig_l1, function(x) any(x %in% myeloid_set))]
gene_prof[, has_lymphoid := sapply(sig_l1, function(x) any(x %in% lymphoid_set))]
gene_prof[, has_tcell_only := sapply(sig_l1, function(x) {
  l1 <- setdiff(x, ORACLE_PBMC)
  length(l1) > 0 && all(l1 %in% ORACLE_TCELL)
})]
gene_prof[, pbmc_only := sapply(sig_l1, function(x) length(x) == 1 && x[1] == ORACLE_PBMC)]

# Stratified sampling
sample_n <- function(dt, n) {
  if (nrow(dt) == 0) return(dt[0])
  dt[sample(.N, min(n, .N))]
}

gene_samples <- rbind(
  sample_n(gene_prof[has_myeloid == TRUE & has_lymphoid == TRUE], 3),     # cross-lineage
  sample_n(gene_prof[has_myeloid != has_lymphoid & n_sig_l1 > 1], 3),     # lineage-specific
  sample_n(gene_prof[has_tcell_only == TRUE & n_sig_l1 > 1], 3),          # T-cell only
  sample_n(gene_prof[n_sig_l1 == 1 & pbmc_only == FALSE], 3),             # single CT
  sample_n(gene_prof[pbmc_only == TRUE], 2),                               # PBMC only
  sample_n(gene_prof[n_sig_l1 == 0 & pbmc_only == FALSE], 3),             # no sig
  sample_n(gene_prof[n_tested >= 30], 3)                                   # high CT count
)
gene_samples <- unique(gene_samples, by = "gene_id")
cat("  Sampled", nrow(gene_samples), "genes\n")

# Load LFSR for sampled genes
eqtl_lfsr <- fread(eqtl_lfsr_file)

# Build fixtures
gene_fixtures <- lapply(seq_len(nrow(gene_samples)), function(i) {
  gid <- gene_samples$gene_id[i]

  # ACAT q-values
  acat <- gene_acat[gene_id == gid, .(cell_type, q_value)]

  # LFSR: find rows matching this gene
  lfsr_rows <- eqtl_lfsr[startsWith(id, paste0(gid, ":"))]
  lfsr_data <- NULL
  if (nrow(lfsr_rows) > 0) {
    # Take the row with the most significant variant (min LFSR across CTs)
    ct_cols <- intersect(names(lfsr_rows), l1_cts)
    if (length(ct_cols) > 0) {
      min_lfsr <- lfsr_rows[, lapply(.SD, min, na.rm = TRUE), .SDcols = ct_cols]
      lfsr_data <- data.table(
        cell_type = ct_cols,
        lfsr = as.numeric(min_lfsr[1, ])
      )
    }
  }

  # Oracle derivation
  sig_cts <- acat[q_value < sig_threshold]$cell_type
  sig_l1 <- unique(oracle_map_to_l1(sig_cts))

  lfsr_named <- NULL
  if (!is.null(lfsr_data)) {
    lfsr_named <- setNames(lfsr_data$lfsr, lfsr_data$cell_type)
  }

  oracle_result <- oracle_categorize_specificity(
    sig_l1_cts = sig_l1,
    lfsr_values = lfsr_named,
    tested_l1_cts = unique(oracle_map_to_l1(acat$cell_type))
  )

  list(
    gene_id = gid,
    acat_qvalues = acat,
    lfsr_values = lfsr_data,
    oracle_trace = list(
      significant_cts_raw = sig_cts,
      significant_l1 = sig_l1,
      n_sig_l1 = length(setdiff(sig_l1, ORACLE_PBMC))
    ),
    expected_category = oracle_result
  )
})
names(gene_fixtures) <- sapply(gene_fixtures, `[[`, "gene_id")

# ============================================================================
# PEAK FIXTURES (N=20)
# ============================================================================
cat("=== Generating peak fixtures ===\n")

peak_acat_list <- list()
for (ct in l1_cts) {
  f <- construct_path(caqtl_acat_pattern, cell_type = ct)
  if (file.exists(f)) {
    dt <- fread(f, select = c(caqtl_acat_id_col, caqtl_acat_q_col))
    setnames(dt, c(caqtl_acat_id_col, caqtl_acat_q_col), c("peak_id", "q_value"))
    dt[, cell_type := ct]
    peak_acat_list[[ct]] <- dt
  }
}
peak_acat <- rbindlist(peak_acat_list)

peak_sig <- peak_acat[q_value < sig_threshold, .(sig_cts = list(cell_type)), by = peak_id]
peak_all <- peak_acat[, .(tested_cts = list(unique(cell_type)), n_tested = uniqueN(cell_type)), by = peak_id]
peak_prof <- merge(peak_all, peak_sig, by = "peak_id", all.x = TRUE)
peak_prof[is.na(sig_cts), sig_cts := list(list(character(0)))]

peak_prof[, sig_l1 := lapply(sig_cts, function(x) unique(oracle_map_to_l1(unlist(x))))]
peak_prof[, n_sig_l1 := sapply(sig_l1, function(x) length(setdiff(x, ORACLE_PBMC)))]
peak_prof[, has_myeloid := sapply(sig_l1, function(x) any(x %in% myeloid_set))]
peak_prof[, has_lymphoid := sapply(sig_l1, function(x) any(x %in% lymphoid_set))]
peak_prof[, has_tcell_only := sapply(sig_l1, function(x) {
  l1 <- setdiff(x, ORACLE_PBMC)
  length(l1) > 0 && all(l1 %in% ORACLE_TCELL)
})]
peak_prof[, pbmc_only := sapply(sig_l1, function(x) length(x) == 1 && x[1] == ORACLE_PBMC)]

peak_samples <- rbind(
  sample_n(peak_prof[has_myeloid == TRUE & has_lymphoid == TRUE], 3),
  sample_n(peak_prof[has_myeloid != has_lymphoid & n_sig_l1 > 1], 3),
  sample_n(peak_prof[has_tcell_only == TRUE & n_sig_l1 > 1], 3),
  sample_n(peak_prof[n_sig_l1 == 1 & pbmc_only == FALSE], 3),
  sample_n(peak_prof[pbmc_only == TRUE], 2),
  sample_n(peak_prof[n_sig_l1 == 0 & pbmc_only == FALSE], 3),
  sample_n(peak_prof[n_tested >= 8], 3)
)
peak_samples <- unique(peak_samples, by = "peak_id")
cat("  Sampled", nrow(peak_samples), "peaks\n")

caqtl_lfsr <- fread(caqtl_lfsr_file)

peak_fixtures <- lapply(seq_len(nrow(peak_samples)), function(i) {
  pid <- peak_samples$peak_id[i]

  acat <- peak_acat[peak_id == pid, .(cell_type, q_value)]

  lfsr_rows <- caqtl_lfsr[startsWith(id, paste0(pid, ":"))]
  lfsr_data <- NULL
  if (nrow(lfsr_rows) > 0) {
    ct_cols <- intersect(names(lfsr_rows), l1_cts)
    if (length(ct_cols) > 0) {
      min_lfsr <- lfsr_rows[, lapply(.SD, min, na.rm = TRUE), .SDcols = ct_cols]
      lfsr_data <- data.table(cell_type = ct_cols, lfsr = as.numeric(min_lfsr[1, ]))
    }
  }

  sig_cts <- acat[q_value < sig_threshold]$cell_type
  sig_l1 <- unique(oracle_map_to_l1(sig_cts))

  lfsr_named <- NULL
  if (!is.null(lfsr_data)) lfsr_named <- setNames(lfsr_data$lfsr, lfsr_data$cell_type)

  oracle_result <- oracle_categorize_specificity(
    sig_l1_cts = sig_l1,
    lfsr_values = lfsr_named,
    tested_l1_cts = unique(oracle_map_to_l1(acat$cell_type))
  )

  list(
    peak_id = pid,
    acat_qvalues = acat,
    lfsr_values = lfsr_data,
    oracle_trace = list(
      significant_cts_raw = sig_cts,
      significant_l1 = sig_l1,
      n_sig_l1 = length(setdiff(sig_l1, ORACLE_PBMC))
    ),
    expected_category = oracle_result
  )
})
names(peak_fixtures) <- sapply(peak_fixtures, `[[`, "peak_id")

# ============================================================================
# VARIANT FIXTURES (N=30)
# ============================================================================
cat("=== Generating variant fixtures ===\n")

# Load production output for stratified sampling by mechanism
# Uses output_dir from config, or a command-line override
prod_output_dir <- config$output_dir %||% "results"
prod_variant_file <- file.path(prod_output_dir, "variant_categorization.tsv.gz")
if (!file.exists(prod_variant_file)) {
  # Fallback: try relative to config file location
  prod_variant_file <- file.path(dirname(config_file), prod_output_dir, "variant_categorization.tsv.gz")
}
if (!file.exists(prod_variant_file)) {
  stop("Cannot find variant_categorization.tsv.gz for sampling. Run production first or set output_dir in config.")
}
prod_variants <- fread(prod_variant_file)

# Sample by mechanism category (using output for ID discovery only, not for expected values)
# Old "Full Cascade" split into three mechanisms: Local, Positional, Distal Cascade.
variant_samples <- rbind(
  sample_n(prod_variants[qtl_mechanism_category == "Local Cascade"], 2),
  sample_n(prod_variants[qtl_mechanism_category == "Positional Cascade"], 2),
  sample_n(prod_variants[qtl_mechanism_category == "Distal Cascade"], 2),
  sample_n(prod_variants[qtl_mechanism_category == "caQTL + eQTL (No Link)"], 4),
  sample_n(prod_variants[qtl_mechanism_category == "Only caQTL (With Link)"], 3),
  sample_n(prod_variants[qtl_mechanism_category == "Only caQTL (No Link)"], 4),
  sample_n(prod_variants[qtl_mechanism_category == "Only eQTL"], 4),
  # Multi-CT diverse: different specificity categories
  sample_n(prod_variants[cell_type_specificity == "Cross-lineage shared"], 3),
  sample_n(prod_variants[cell_type_specificity == "Likely shared but underpowered"], 3),
  sample_n(prod_variants[cell_type_specificity == "Single cell-type" & qtl_mechanism_category == "Local Cascade"], 2),
  sample_n(prod_variants[cell_type_specificity == "T-cell-specific"], 2)
)
variant_samples <- unique(variant_samples, by = "variant_id")
cat("  Sampled", nrow(variant_samples), "variants\n")

# For each sampled variant, extract per-CT SuSiE data and derive expected pattern
# This is the most complex part — we need to trace through raw data
variant_fixtures <- lapply(seq_len(nrow(variant_samples)), function(i) {
  vid <- variant_samples$variant_id[i]
  chr <- sub("_.*", "", vid)

  # Store the production output for comparison (not as expected — oracle derives expected)
  prod_row <- prod_variants[variant_id == vid]

  # Extract per-CT data — scan ALL cell types (L1 + L2) including PBMC
  # The production pipeline processes all 33 CTs; L2-only or PBMC-only eQTL would be missed
  per_ct <- list()
  for (ct in cell_types) {
    # eQTL SuSiE
    eqtl_f <- construct_path(config$file_patterns$eqtl_susie, cell_type = ct, chr = chr)
    eqtl_pip <- 0
    eqtl_genes <- character(0)
    if (file.exists(eqtl_f) && file.size(eqtl_f) > 0) {
      eqtl <- tryCatch({
        dt <- fread(eqtl_f, select = c("rsid", "region", "prob"))
        if (nrow(dt) > 0) dt else data.table()
      }, error = function(e) data.table())
      if (nrow(eqtl) > 0 && "rsid" %in% names(eqtl)) {
        eqtl_hit <- eqtl[rsid == vid & prob > 0.01]
        if (nrow(eqtl_hit) > 0) {
          eqtl_pip <- max(eqtl_hit$prob)
          eqtl_genes <- unique(eqtl_hit$region)
        }
      }
    }

    # caQTL SuSiE
    caqtl_f <- construct_path(config$file_patterns$caqtl_susie, cell_type = ct, chr = chr)
    caqtl_pip <- 0
    caqtl_peaks <- character(0)
    if (file.exists(caqtl_f) && file.size(caqtl_f) > 0) {
      caqtl <- tryCatch({
        dt <- fread(caqtl_f, select = c("rsid", "region", "prob"))
        if (nrow(dt) > 0) dt else data.table()
      }, error = function(e) data.table())
      if (nrow(caqtl) > 0 && "rsid" %in% names(caqtl)) {
        caqtl_hit <- caqtl[rsid == vid & prob > 0.01]
        if (nrow(caqtl_hit) > 0) {
          caqtl_pip <- max(caqtl_hit$prob)
          caqtl_peaks <- unique(caqtl_hit$region)
        }
      }
    }

    per_ct[[ct]] <- list(
      eqtl_pip = eqtl_pip,
      eqtl_genes = eqtl_genes,
      caqtl_pip = caqtl_pip,
      caqtl_peaks = caqtl_peaks
    )
  }

  list(
    variant_id = vid,
    chromosome = chr,
    per_celltype = per_ct,
    production_output = list(
      qtl_mechanism_category = prod_row$qtl_mechanism_category,
      cell_type_specificity = prod_row$cell_type_specificity,
      qtl_pattern_number = prod_row$qtl_pattern_number,
      peak_gene_link = prod_row$peak_gene_link,
      eqtl = prod_row$eqtl
    )
  )
})
names(variant_fixtures) <- sapply(variant_fixtures, `[[`, "variant_id")

# ============================================================================
# Save fixtures
# ============================================================================
cat("\n=== Saving fixtures ===\n")

saveRDS(gene_fixtures, file.path(output_dir, "validation_genes.rds"))
cat("  Gene fixtures:", length(gene_fixtures), "genes,",
    format(file.size(file.path(output_dir, "validation_genes.rds")), big.mark = ","), "bytes\n")

saveRDS(peak_fixtures, file.path(output_dir, "validation_peaks.rds"))
cat("  Peak fixtures:", length(peak_fixtures), "peaks,",
    format(file.size(file.path(output_dir, "validation_peaks.rds")), big.mark = ","), "bytes\n")

saveRDS(variant_fixtures, file.path(output_dir, "validation_variants.rds"))
cat("  Variant fixtures:", length(variant_fixtures), "variants,",
    format(file.size(file.path(output_dir, "validation_variants.rds")), big.mark = ","), "bytes\n")

cat("\nDone!\n")

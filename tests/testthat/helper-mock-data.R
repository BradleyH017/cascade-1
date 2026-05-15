# Helper functions for creating mock data for CASCADE tests
library(data.table)

#' Create comprehensive mock data files matching mydata_test_config.json structure
#'
#' @param test_dir Directory to create mock files in
#' @param cell_types Vector of cell type names
#' @param chromosomes Vector of chromosome names
#' @param n_genes Number of genes to simulate
#' @param n_peaks Number of peaks to simulate
#' @param n_variants Number of variants to simulate
#' @return List of file patterns pointing to created mock files
create_comprehensive_mock_data <- function(test_dir,
                                           cell_types = NULL,
                                           chromosomes = NULL,
                                           n_genes = 100,
                                           n_peaks = 50,
                                           n_variants = 200) {
  # Use default cell types if not provided (subset from mydata_test_config.json)
  if (is.null(cell_types)) {
    cell_types <- c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l1.CD4_T",
      "predicted.celltype.l1.CD8_T",
      "predicted.celltype.l1.B",
      "predicted.celltype.l1.NK"
    )
  }

  # Use default chromosomes if not provided
  if (is.null(chromosomes)) {
    chromosomes <- c("chr7", "chr22")
  }

  # Create directory structure
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)

  # Generate gene and peak IDs
  gene_ids <- paste0("ENSG", sprintf("%011d", 1:n_genes))
  # Use underscore instead of colon to avoid parsing issues
  # Generate start positions and ensure end > start
  start_positions <- sample(1000000:50000000, n_peaks)
  peak_ids <- paste0(
    "chr", sample(c(1:22, "X", "Y"), n_peaks, replace = TRUE),
    "_", start_positions, "_",
    start_positions + 1000 # Ensure end is always 1000 bp after start
  )
  # Create variant IDs in the format expected by the package: chr_pos_ref_alt
  variant_chrs <- sample(c(paste0("chr", 1:22), "chrX", "chrY"), n_variants, replace = TRUE)
  variant_positions <- sample(1000000:50000000, n_variants)
  variant_refs <- sample(c("A", "T", "C", "G"), n_variants, replace = TRUE)
  variant_alts <- sample(c("A", "T", "C", "G"), n_variants, replace = TRUE)
  # Ensure ref != alt
  for (i in 1:length(variant_alts)) {
    if (variant_refs[i] == variant_alts[i]) {
      variant_alts[i] <- sample(setdiff(c("A", "T", "C", "G"), variant_refs[i]), 1)
    }
  }
  variant_ids <- paste0(variant_chrs, "_", variant_positions, "_", variant_refs, "_", variant_alts)

  # 1. Create eQTL ACAT files
  for (ct in cell_types) {
    acat_data <- data.table(
      phenotype_id = gene_ids,
      variant_id = sample(variant_ids, n_genes, replace = TRUE),
      pval = runif(n_genes, 1e-10, 1),
      ACAT_p = runif(n_genes, 1e-10, 1),
      ACAT_q = p.adjust(runif(n_genes, 0, 1), method = "BH"),
      beta = rnorm(n_genes, 0, 0.5),
      se = runif(n_genes, 0.01, 0.1)
    )

    acat_file <- file.path(test_dir, gsub(
      "\\{CELL_TYPE\\}", ct,
      "integrated_gex_batch1_5.fgid.qc.{CELL_TYPE}.mean.inv.SAIGE.acat.txt.gz"
    ))
    dir.create(dirname(acat_file), recursive = TRUE, showWarnings = FALSE)
    fwrite(acat_data, acat_file, sep = "\t", compress = "gzip")
  }

  # 2. Create caQTL ACAT files
  for (ct in cell_types) {
    caqtl_data <- data.table(
      phenotype_id = peak_ids,
      variant_id = sample(variant_ids, n_peaks, replace = TRUE),
      pval = runif(n_peaks, 1e-10, 1),
      qval = p.adjust(runif(n_peaks, 0, 1), method = "BH"),
      beta = rnorm(n_peaks, 0, 0.3),
      se = runif(n_peaks, 0.01, 0.08)
    )

    caqtl_file <- file.path(test_dir, gsub(
      "\\{CELL_TYPE\\}", ct,
      "integrated_atac_batch1_5.fgid.{CELL_TYPE}.sum.inv.cis_qtl_acat.tsv.gz"
    ))
    dir.create(dirname(caqtl_file), recursive = TRUE, showWarnings = FALSE)
    fwrite(caqtl_data, caqtl_file, sep = "\t", compress = "gzip")
  }

  # 3. Create eQTL SuSiE files
  susie_dir <- file.path(test_dir, "eqtl_susie")
  dir.create(susie_dir, recursive = TRUE, showWarnings = FALSE)

  for (ct in cell_types) {
    for (chr in chromosomes) {
      # Create credible sets data
      n_cs_variants <- min(n_variants, 50)
      # Use the same variant IDs for consistency
      susie_data <- data.table(
        rsid = sample(variant_ids, n_cs_variants),
        region = sample(gene_ids, n_cs_variants, replace = TRUE),
        prob = runif(n_cs_variants, 0, 1),
        cs = sample(1:5, n_cs_variants, replace = TRUE),
        chromosome = chr,
        position = sample(1000000:50000000, n_cs_variants),
        ref = sample(c("A", "T", "C", "G"), n_cs_variants, replace = TRUE),
        alt = sample(c("A", "T", "C", "G"), n_cs_variants, replace = TRUE)
      )
      # Normalize PIP values within credible sets
      susie_data[, prob := prob / sum(prob), by = .(region, cs)]

      susie_file <- file.path(susie_dir, paste0(
        "integrated_gex_batch1_5.fgid.qc.", ct, ".mean.inv.SAIGE.",
        chr, ".SUSIE.in_cs.snp.bgz"
      ))
      fwrite(susie_data, susie_file, sep = "\t", compress = "gzip")
    }
  }

  # 4. Create caQTL SuSiE files
  caqtl_susie_dir <- file.path(test_dir, "caqtl_susie")
  dir.create(caqtl_susie_dir, recursive = TRUE, showWarnings = FALSE)

  for (ct in cell_types) {
    for (chr in chromosomes) {
      n_cs_variants <- min(n_variants, 40)
      caqtl_susie_data <- data.table(
        rsid = sample(variant_ids, n_cs_variants), # Using same format as eQTL
        region = sample(peak_ids, n_cs_variants, replace = TRUE),
        prob = runif(n_cs_variants, 0, 1),
        cs = sample(1:3, n_cs_variants, replace = TRUE),
        chromosome = chr,
        position = sample(1000000:50000000, n_cs_variants)
      )
      caqtl_susie_data[, prob := prob / sum(prob), by = .(region, cs)]

      caqtl_susie_file <- file.path(caqtl_susie_dir, paste0(
        "integrated_atac_batch1_5.fgid.", ct, ".sum.inv.",
        chr, ".SUSIE.in_cs.snp.bgz"
      ))
      fwrite(caqtl_susie_data, caqtl_susie_file, sep = "\t", compress = "gzip")
    }
  }

  # 5. Create peak-gene links files
  for (ct in cell_types) {
    peak_gene_data <- data.table(
      peak_id = sample(peak_ids, min(n_peaks * 2, length(peak_ids) * 2), replace = TRUE),
      gene_id = sample(gene_ids, min(n_peaks * 2, length(gene_ids) * 2), replace = TRUE),
      Score = runif(min(n_peaks * 2, length(peak_ids) * 2), 0, 10),
      FDR = p.adjust(runif(min(n_peaks * 2, length(peak_ids) * 2), 0, 1), method = "BH")
    )

    pg_file <- file.path(test_dir, paste0("open4gene.", ct, ".results.sig.tsv.gz"))
    fwrite(peak_gene_data, pg_file, sep = "\t", compress = "gzip")
  }

  # 6. Create peak BED files
  for (ct in cell_types) {
    # Parse peak IDs with underscore format
    peak_parts <- strsplit(peak_ids, "_")
    bed_data <- data.table(
      chr = sapply(peak_parts, `[`, 1),
      start = as.integer(sapply(peak_parts, `[`, 2)),
      end = as.integer(sapply(peak_parts, `[`, 3)),
      peak_id = peak_ids,
      score = rep(1000, length(peak_ids)), # Standard BED format requires score
      strand = rep("*", length(peak_ids)) # Standard BED format requires strand
    )

    bed_file <- file.path(test_dir, paste0(
      "integrated_atac_batch1_5.fgid.", ct,
      ".sum.inv.cis_qtl_pairs.peak.txt"
    ))
    fwrite(bed_data, bed_file, sep = "\t", col.names = FALSE)
  }

  # 7. Create LFSR files (critical for many tests)
  # eQTL LFSR
  eqtl_lfsr_data <- data.table()
  for (i in 1:min(n_genes * 5, 500)) {
    gene_id <- sample(gene_ids, 1)
    variant_id <- sample(variant_ids, 1)
    row_data <- data.table(
      id = paste0(gene_id, ":", variant_id), # Required format: feature_id:variant_id
      variant = variant_id,
      gene = gene_id
    )
    # Add LFSR values for each cell type
    for (ct in cell_types) {
      # Create realistic LFSR values (mostly high, some low for significance)
      if (runif(1) < 0.2) { # 20% significant
        row_data[[paste0("lfsr_", ct)]] <- runif(1, 0, 0.05)
      } else {
        row_data[[paste0("lfsr_", ct)]] <- runif(1, 0.5, 1)
      }
    }
    eqtl_lfsr_data <- rbind(eqtl_lfsr_data, row_data)
  }

  eqtl_lfsr_file <- file.path(test_dir, "integrated_gex_batch1_5.mashr.lfsr.tsv.gz")
  fwrite(eqtl_lfsr_data, eqtl_lfsr_file, sep = "\t", compress = "gzip")

  # caQTL LFSR
  caqtl_lfsr_data <- data.table()
  for (i in 1:min(n_peaks * 5, 300)) {
    peak_id <- sample(peak_ids, 1)
    variant_id <- sample(variant_ids, 1)
    row_data <- data.table(
      id = paste0(peak_id, ":", variant_id), # Required format: feature_id:variant_id
      variant = variant_id,
      peak = peak_id
    )
    for (ct in cell_types) {
      if (runif(1) < 0.15) { # 15% significant
        row_data[[paste0("lfsr_", ct)]] <- runif(1, 0, 0.05)
      } else {
        row_data[[paste0("lfsr_", ct)]] <- runif(1, 0.5, 1)
      }
    }
    caqtl_lfsr_data <- rbind(caqtl_lfsr_data, row_data)
  }

  caqtl_lfsr_file <- file.path(test_dir, "integrated_atac_batch1_5.mashr.lfsr.tsv.gz")
  fwrite(caqtl_lfsr_data, caqtl_lfsr_file, sep = "\t", compress = "gzip")

  # 8. Create meta data files (required for Cochran's Q analysis)
  # eQTL meta data
  n_meta_eqtl <- min(n_genes * 3, 150)
  eqtl_meta_data <- data.table(
    variant = sample(variant_ids, n_meta_eqtl, replace = TRUE),
    phenotype = sample(gene_ids, n_meta_eqtl, replace = TRUE),
    meta_nlog10p_het = runif(n_meta_eqtl, 0, 20),
    direction = sample(c("+", "-", "+-", "-+"), n_meta_eqtl, replace = TRUE),
    max_pip = runif(n_meta_eqtl, 0, 1),
    max_chisq = runif(n_meta_eqtl, 0, 50)
  )
  eqtl_meta_file <- file.path(test_dir, "eqtl_meta_data.tsv.gz")
  fwrite(eqtl_meta_data, eqtl_meta_file, sep = "\t", compress = "gzip")

  # caQTL meta data
  n_meta_caqtl <- min(n_peaks * 3, 90)
  caqtl_meta_data <- data.table(
    variant = sample(variant_ids, n_meta_caqtl, replace = TRUE),
    phenotype = sample(peak_ids, n_meta_caqtl, replace = TRUE),
    meta_nlog10p_het = runif(n_meta_caqtl, 0, 20),
    direction = sample(c("+", "-", "+-", "-+"), n_meta_caqtl, replace = TRUE),
    max_pip = runif(n_meta_caqtl, 0, 1),
    max_chisq = runif(n_meta_caqtl, 0, 50)
  )
  caqtl_meta_file <- file.path(test_dir, "caqtl_meta_data.tsv.gz")
  fwrite(caqtl_meta_data, caqtl_meta_file, sep = "\t", compress = "gzip")

  # 9. Create mashr model RDS files (optional but useful)
  # Create simplified mock mashr models
  mock_mashr_eqtl <- list(
    fitted_g = list(
      Ulist = list(
        tFlash_1 = diag(length(cell_types)),
        tFlash_2 = matrix(0.5, length(cell_types), length(cell_types))
      ),
      grid = c(0.5, 1, 2, 4, 8),
      pi = rep(0.2, 10)
    ),
    result = list(
      PosteriorMean = matrix(rnorm(100 * length(cell_types)), 100, length(cell_types)),
      PosteriorSD = matrix(runif(100 * length(cell_types), 0.1, 1), 100, length(cell_types)),
      lfsr = matrix(runif(100 * length(cell_types), 0, 1), 100, length(cell_types))
    )
  )

  mock_mashr_caqtl <- mock_mashr_eqtl # Similar structure for caQTL

  eqtl_mashr_file <- file.path(test_dir, "integrated_gex_batch1_5.mashr.model.strong.rds")
  saveRDS(mock_mashr_eqtl, eqtl_mashr_file)

  caqtl_mashr_file <- file.path(test_dir, "integrated_atac_batch1_5.mashr.model.strong.rds")
  saveRDS(mock_mashr_caqtl, caqtl_mashr_file)

  # Return file patterns matching mydata_test_config.json structure
  list(
    eqtl_acat = file.path(test_dir, "integrated_gex_batch1_5.fgid.qc.{CELL_TYPE}.mean.inv.SAIGE.acat.txt.gz"),
    caqtl_acat = file.path(test_dir, "integrated_atac_batch1_5.fgid.{CELL_TYPE}.sum.inv.cis_qtl_acat.tsv.gz"),
    eqtl_susie = file.path(test_dir, "eqtl_susie/integrated_gex_batch1_5.fgid.qc.{CELL_TYPE}.mean.inv.SAIGE.{CHR}.SUSIE.in_cs.snp.bgz"),
    caqtl_susie = file.path(test_dir, "caqtl_susie/integrated_atac_batch1_5.fgid.{CELL_TYPE}.sum.inv.{CHR}.SUSIE.in_cs.snp.bgz"),
    peak_gene_links = file.path(test_dir, "open4gene.{CELL_TYPE}.results.sig.tsv.gz"),
    peak_bed = file.path(test_dir, "integrated_atac_batch1_5.fgid.{CELL_TYPE}.sum.inv.cis_qtl_pairs.peak.txt"),
    eqtl_lfsr = eqtl_lfsr_file,
    caqtl_lfsr = caqtl_lfsr_file,
    eqtl_mashr = eqtl_mashr_file,
    caqtl_mashr = caqtl_mashr_file,
    eqtl_meta = eqtl_meta_file,
    caqtl_meta = caqtl_meta_file
  )
}

#' Create a test configuration matching mydata_test_config.json structure
#'
#' @param test_dir Directory for test files
#' @param cell_types Cell types to include
#' @param chromosomes Chromosomes to include
#' @param feature_type Type of features to analyze
#' @return Configuration list
create_test_config <- function(test_dir,
                               cell_types = NULL,
                               chromosomes = NULL,
                               feature_type = "all") {
  # Use defaults if not provided
  if (is.null(cell_types)) {
    cell_types <- c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l1.CD4_T",
      "predicted.celltype.l1.CD8_T",
      "predicted.celltype.l1.B",
      "predicted.celltype.l1.NK"
    )
  }

  if (is.null(chromosomes)) {
    chromosomes <- c("chr7")
  }

  # Create mock data files
  file_patterns <- create_comprehensive_mock_data(
    test_dir = test_dir,
    cell_types = cell_types,
    chromosomes = chromosomes,
    n_genes = 50,
    n_peaks = 30,
    n_variants = 100
  )

  # Build configuration matching mydata_test_config.json structure
  config <- list(
    cell_types = cell_types,
    chromosomes = chromosomes,
    feature_type = feature_type,
    file_patterns = file_patterns,
    parameters = list(
      n_cores = 1,
      pip_threshold = 0.5,
      acat_fdr_threshold = 0.05,
      lfsr_sig_threshold = 0.05,
      lfsr_null_threshold = 0.5,
      use_cochran_q = TRUE,
      cochran_q_threshold = 5e-8,
      peak_gene_link_fdr_threshold = 0.05
    ),
    column_mapping = list(
      eqtl_acat = list(
        feature_id = "phenotype_id",
        q_value = "ACAT_q"
      ),
      caqtl_acat = list(
        feature_id = "phenotype_id",
        q_value = "qval"
      ),
      eqtl_susie = list(
        variant_id = "rsid",
        feature_id = "region",
        pip = "prob",
        chromosome = "chromosome",
        cs_id = "cs",
        beta = NULL,
        se = NULL
      ),
      caqtl_susie = list(
        variant_id = "rsid",
        feature_id = "region",
        pip = "prob",
        chromosome = "chromosome",
        cs_id = "cs",
        beta = NULL,
        se = NULL
      )
    ),
    debug = list(
      enabled = FALSE,
      n_genes = 50,
      n_peaks = 30
    ),
    cache = list(
      enabled = TRUE
    ),
    output_dir = test_dir
  )

  return(config)
}

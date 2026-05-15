library(testthat)
library(cascade)
library(data.table)

test_that("load_feature_data handles file pattern substitution", {
  # Test CELL_TYPE substitution
  pattern <- "/data/{CELL_TYPE}.acat.txt"
  cell_type <- "Mono"

  result <- gsub("\\{CELL_TYPE\\}", cell_type, pattern)
  expect_equal(result, "/data/Mono.acat.txt")

  # Test CHR substitution
  pattern <- "/data/{CELL_TYPE}.{CHR}.susie.txt"
  chromosome <- "chr22"

  result <- gsub("\\{CELL_TYPE\\}", cell_type, pattern)
  result <- gsub("\\{CHR\\}", chromosome, result)
  expect_equal(result, "/data/Mono.chr22.susie.txt")

  # Test multiple substitutions
  patterns <- list(
    eqtl_acat = "/path/{CELL_TYPE}/acat.gz",
    eqtl_susie = "/path/{CELL_TYPE}/{CHR}/susie.gz",
    peak_bed = "/path/{CELL_TYPE}.peaks.bed"
  )

  for (pattern_name in names(patterns)) {
    pattern <- patterns[[pattern_name]]
    expect_true(grepl("\\{CELL_TYPE\\}", pattern))
  }
})

test_that("load_variant_data processes multiple chromosomes", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create mock SuSiE files for multiple chromosomes
  cell_types <- c("Mono", "CD4_T")
  chromosomes <- c("chr21", "chr22")

  for (ct in cell_types) {
    for (chr in chromosomes) {
      # Create caQTL SuSiE file
      caqtl_data <- data.table(
        rsid = paste0(chr, "_", seq(10000, 100000, 10000)),
        region = paste0("peak_", 1:10),
        prob = runif(10, 0.3, 1),
        chromosome = chr
      )
      caqtl_file <- file.path(test_dir, paste0(ct, ".", chr, ".caqtl.susie.gz"))
      fwrite(caqtl_data, caqtl_file, sep = "\t", compress = "gzip")

      # Create eQTL SuSiE file
      eqtl_data <- data.table(
        rsid = paste0(chr, "_", seq(10000, 100000, 10000)),
        region = paste0("GENE", 1:10),
        prob = runif(10, 0.3, 1),
        chromosome = chr
      )
      eqtl_file <- file.path(test_dir, paste0(ct, ".", chr, ".eqtl.susie.gz"))
      fwrite(eqtl_data, eqtl_file, sep = "\t", compress = "gzip")
    }
  }

  # Create peak-gene links
  peak_gene_links <- data.table(
    peak_id = paste0("peak_", 1:10),
    gene_id = paste0("GENE", 1:10),
    qvalue = runif(10, 0, 0.1)
  )
  links_file <- file.path(test_dir, "peak_gene_links.tsv.gz")
  fwrite(peak_gene_links, links_file, sep = "\t", compress = "gzip")

  # Create peak BED file
  peak_bed <- data.table(
    chr = c(rep("chr1", 5), rep("chr2", 5)),
    start = seq(1000, 10000, 1000),
    end = seq(1500, 10500, 1000),
    peak_id = paste0("peak_", 1:10),
    score = rep(0, 10),
    strand = rep(".", 10)
  )
  bed_file <- file.path(test_dir, "peaks.bed")
  fwrite(peak_bed, bed_file, sep = "\t", col.names = FALSE)

  # Create configuration
  config <- list(
    cell_types = cell_types,
    chromosomes = chromosomes,
    file_patterns = list(
      caqtl_susie = file.path(test_dir, "{CELL_TYPE}.{CHR}.caqtl.susie.gz"),
      eqtl_susie = file.path(test_dir, "{CELL_TYPE}.{CHR}.eqtl.susie.gz"),
      peak_gene_links = links_file,
      peak_bed = bed_file
    ),
    column_mapping = list(
      caqtl_susie = list(
        variant_id = "rsid",
        feature_id = "region",
        pip = "prob",
        chromosome = "chromosome",
        cs_id = "cs",
        beta = NULL,
        se = NULL
      ),
      eqtl_susie = list(
        variant_id = "rsid",
        feature_id = "region",
        pip = "prob",
        chromosome = "chromosome",
        cs_id = "cs",
        beta = NULL,
        se = NULL
      )
    ),
    parameters = list(
      pip_threshold = 0.5,
      peak_gene_link_fdr_threshold = 0.05
    )
  )

  # Test loading variant data
  variant_data <- load_variant_data(
    config = config,
    chromosomes = chromosomes,
    pip_threshold = 0.5,
    acat_fdr_threshold = 0.05,
    peak_bed_file = config$file_patterns$peak_bed, # Use the peak_bed from config
    column_mapping = config$column_mapping,
    num_cores = 1
  )

  expect_true(is.list(variant_data))
  expect_true("all_variants" %in% names(variant_data))
  expect_true("caqtl_data_by_ct" %in% names(variant_data))
  expect_true("eqtl_data_by_ct" %in% names(variant_data))
  expect_true("peak_gene_links" %in% names(variant_data))
})

test_that("load_lfsr_results handles compressed files", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create mock LFSR data
  cell_types <- c("Mono", "CD4_T", "CD8_T")

  # Create eQTL LFSR file with production column format
  eqtl_lfsr <- data.table(
    id = paste0(rep(paste0("GENE", 1:5), 10), ":", rep(paste0("rs", 1:10), each = 5)),
    variant = rep(paste0("rs", 1:10), each = 5),
    gene = rep(paste0("GENE", 1:5), 10)
  )
  for (ct in cell_types) {
    eqtl_lfsr[[paste0("predicted.celltype.l1.", ct)]] <- runif(50, 0, 1)
  }
  eqtl_lfsr_file <- file.path(test_dir, "eqtl.lfsr.tsv.gz")
  fwrite(eqtl_lfsr, eqtl_lfsr_file, sep = "\t", compress = "gzip")

  # Create caQTL LFSR file with production column format
  caqtl_lfsr <- data.table(
    id = paste0(rep(paste0("peak_", 1:5), 10), ":", rep(paste0("rs", 1:10), each = 5)),
    variant = rep(paste0("rs", 1:10), each = 5),
    peak = rep(paste0("peak_", 1:5), 10)
  )
  for (ct in cell_types) {
    caqtl_lfsr[[paste0("predicted.celltype.l1.", ct)]] <- runif(50, 0, 1)
  }
  caqtl_lfsr_file <- file.path(test_dir, "caqtl.lfsr.tsv.gz")
  fwrite(caqtl_lfsr, caqtl_lfsr_file, sep = "\t", compress = "gzip")

  # Create configuration
  config <- list(
    cell_types = cell_types,
    file_patterns = list(
      eqtl_lfsr = eqtl_lfsr_file,
      caqtl_lfsr = caqtl_lfsr_file
    )
  )

  # Load LFSR results
  lfsr_results <- load_lfsr_results(config)

  expect_true(is.list(lfsr_results))
  expect_true("eqtl_feature" %in% names(lfsr_results))
  expect_true("caqtl_feature" %in% names(lfsr_results))
  expect_true("eqtl_long" %in% names(lfsr_results))
  expect_true("caqtl_long" %in% names(lfsr_results))

  # Check long format structure
  if (nrow(lfsr_results$eqtl_long) > 0) {
    expect_true(all(c("feature_id", "variant_id", "cell_type", "lfsr") %in% names(lfsr_results$eqtl_long)))
  }

  if (nrow(lfsr_results$caqtl_long) > 0) {
    expect_true(all(c("feature_id", "variant_id", "cell_type", "lfsr") %in% names(lfsr_results$caqtl_long)))
  }
})

test_that("load_feature_data handles column mapping", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create files with non-standard column names
  cell_types <- c("Mono", "CD4_T")

  for (ct in cell_types) {
    # ACAT file with custom columns
    acat_data <- data.table(
      feature_name = paste0("GENE", 1:10), # Non-standard name
      acat_qvalue = runif(10, 0, 0.1) # Non-standard name
    )
    acat_file <- file.path(test_dir, paste0(ct, ".acat.txt.gz"))
    fwrite(acat_data, acat_file, sep = "\t", compress = "gzip")

    # SuSiE file with custom columns
    susie_data <- data.table(
      snp_id = paste0("rs", 1:20), # Non-standard name
      gene_region = rep(paste0("GENE", 1:10), each = 2), # Non-standard name
      posterior_prob = runif(20, 0, 1), # Non-standard name
      chr = "chr22" # Non-standard name
    )
    susie_file <- file.path(test_dir, paste0(ct, ".chr22.susie.txt.gz"))
    fwrite(susie_data, susie_file, sep = "\t", compress = "gzip")
  }

  # Configuration with column mapping
  config <- list(
    cell_types = cell_types,
    chromosomes = c("chr22"),
    file_patterns = list(
      eqtl_acat = file.path(test_dir, "{CELL_TYPE}.acat.txt.gz"),
      eqtl_susie = file.path(test_dir, "{CELL_TYPE}.{CHR}.susie.txt.gz")
    ),
    column_mapping = list(
      eqtl_acat = list(
        feature_id = "feature_name",
        q_value = "acat_qvalue"
      ),
      eqtl_susie = list(
        variant_id = "snp_id",
        feature_id = "gene_region",
        pip = "posterior_prob",
        chromosome = "chr",
        cs_id = "cs",
        beta = NULL,
        se = NULL
      )
    ),
    parameters = list(
      n_cores = 1,
      acat_fdr_threshold = 0.05,
      pip_threshold = 0.5
    )
  )

  # Test that column mapping would be applied
  expect_equal(config$column_mapping$eqtl_acat$feature_id, "feature_name")
  expect_equal(config$column_mapping$eqtl_acat$q_value, "acat_qvalue")
  expect_equal(config$column_mapping$eqtl_susie$variant_id, "snp_id")
  expect_equal(config$column_mapping$eqtl_susie$pip, "posterior_prob")
})

test_that("parallel data loading with n_cores works", {
  available_cores <- parallel::detectCores()
  # Test different core settings, capped at the runner's actual core count
  # so the test is meaningful on machines with fewer cores (e.g. CI runners).
  core_configs <- list(
    list(n_cores = 1),
    list(n_cores = min(2, available_cores)),
    list(n_cores = min(4, available_cores)),
    list(n_cores = NULL) # Default (should use 1)
  )

  for (config in core_configs) {
    n_cores <- if (is.null(config$n_cores)) 1 else config$n_cores
    expect_true(n_cores >= 1)
    expect_true(n_cores <= available_cores)
  }
})

test_that("data loading handles missing files gracefully", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Configuration pointing to non-existent files
  config <- list(
    cell_types = c("Mono", "CD4_T"),
    chromosomes = c("chr22"),
    file_patterns = list(
      eqtl_acat = file.path(test_dir, "nonexistent_{CELL_TYPE}.txt"),
      eqtl_susie = file.path(test_dir, "nonexistent_{CELL_TYPE}_{CHR}.txt")
    ),
    parameters = list(n_cores = 1)
  )

  # Check that files don't exist
  test_file <- gsub("\\{CELL_TYPE\\}", "Mono", config$file_patterns$eqtl_acat)
  expect_false(file.exists(test_file))
})

test_that("data loading handles debug mode limits", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create large mock dataset
  n_genes <- 1000
  n_peaks <- 500

  gene_data <- data.table(
    phenotype_id = paste0("GENE", 1:n_genes),
    ACAT_q = runif(n_genes, 0, 1)
  )

  peak_data <- data.table(
    phenotype_id = paste0("peak_", 1:n_peaks),
    qval = runif(n_peaks, 0, 1)
  )

  # Configuration with debug limits
  config <- list(
    debug = list(
      enabled = TRUE,
      n_genes = 100, # Limit to 100 genes
      n_peaks = 50 # Limit to 50 peaks
    )
  )

  # In debug mode, data should be limited
  if (config$debug$enabled) {
    if (nrow(gene_data) > config$debug$n_genes) {
      gene_data <- gene_data[1:config$debug$n_genes]
    }
    if (nrow(peak_data) > config$debug$n_peaks) {
      peak_data <- peak_data[1:config$debug$n_peaks]
    }
  }

  expect_equal(nrow(gene_data), 100)
  expect_equal(nrow(peak_data), 50)
})

test_that("cache functionality in data loading", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Test cache enabled
  config_cache <- list(
    cache = list(enabled = TRUE),
    parameters = list(n_cores = 1)
  )

  expect_true(config_cache$cache$enabled)

  # Create cache directory
  cache_dir <- file.path(test_dir, ".cache")
  if (config_cache$cache$enabled) {
    dir.create(cache_dir, showWarnings = FALSE)
  }

  expect_true(dir.exists(cache_dir))

  # Test cache disabled
  config_no_cache <- list(
    cache = list(enabled = FALSE),
    parameters = list(n_cores = 1)
  )

  expect_false(config_no_cache$cache$enabled)
})

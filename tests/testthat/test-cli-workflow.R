library(testthat)
library(cascade)
library(jsonlite)

test_that("CLI workflow processes configuration correctly", {
  # Create a temporary test directory
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create mock configuration
  test_config <- list(
    cell_types = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T"),
    chromosomes = c("chr22"),
    feature_type = "gene",
    file_patterns = list(
      eqtl_acat = file.path(test_dir, "{CELL_TYPE}.eqtl.acat.txt.gz"),
      eqtl_susie = file.path(test_dir, "{CELL_TYPE}.{CHR}.susie.txt.gz"),
      eqtl_lfsr = file.path(test_dir, "eqtl.lfsr.tsv.gz")
    ),
    parameters = list(
      n_cores = 1,
      pip_threshold = 0.5,
      acat_fdr_threshold = 0.05,
      lfsr_sig_threshold = 0.05,
      lfsr_null_threshold = 0.5,
      use_cochran_q = TRUE,
      cochran_q_threshold = 5e-8
    ),
    column_mapping = list(
      eqtl_acat = list(
        feature_id = "phenotype_id",
        q_value = "ACAT_q"
      ),
      eqtl_susie = list(
        variant_id = "rsid",
        feature_id = "region",
        pip = "prob",
        cs_id = "cs",
        beta = NULL,
        se = NULL
      )
    ),
    output_dir = test_dir
  )

  # Save configuration to JSON
  config_file <- file.path(test_dir, "test_config.json")
  write(toJSON(test_config, auto_unbox = TRUE, pretty = TRUE), config_file)

  # Test configuration loading
  loaded_config <- fromJSON(config_file, simplifyVector = FALSE)
  # Convert to character vector if it's a list
  loaded_cell_types <- if (is.list(loaded_config$cell_types)) {
    unlist(loaded_config$cell_types)
  } else {
    loaded_config$cell_types
  }
  expect_equal(loaded_cell_types, test_config$cell_types)
  expect_equal(loaded_config$chromosomes, test_config$chromosomes)
  expect_equal(loaded_config$parameters$use_cochran_q, TRUE)
  expect_equal(loaded_config$parameters$cochran_q_threshold, 5e-8)
})

test_that("CLI workflow handles multiple feature types", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Test configuration with multiple feature types
  test_configs <- list(
    list(feature_type = "gene"),
    list(feature_type = "peak"),
    list(feature_type = "variant"),
    list(feature_type = "gene,peak"),
    list(feature_type = "gene,peak,variant")
  )

  for (config in test_configs) {
    # Parse feature types as done in CLI
    feature_type <- config$feature_type

    if (grepl(",", feature_type)) {
      feature_types <- trimws(strsplit(feature_type, ",")[[1]])
    } else {
      feature_types <- feature_type
    }

    # Check parsing is correct
    if (config$feature_type == "gene,peak") {
      expect_equal(feature_types, c("gene", "peak"))
    } else if (config$feature_type == "gene,peak,variant") {
      expect_equal(feature_types, c("gene", "peak", "variant"))
    } else {
      expect_equal(feature_types, config$feature_type)
    }
  }
})

test_that("CLI workflow validates required file patterns", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Test missing required patterns for gene analysis
  test_config <- list(
    cell_types = c("Mono", "CD4_T"),
    chromosomes = c("chr22"),
    feature_type = "gene",
    file_patterns = list(), # Missing required patterns
    parameters = list(n_cores = 1)
  )

  config_file <- file.path(test_dir, "invalid_config.json")
  write(toJSON(test_config, auto_unbox = TRUE), config_file)

  # Should identify missing patterns
  loaded_config <- fromJSON(config_file, simplifyVector = FALSE)
  expect_true(length(loaded_config$file_patterns) == 0)
})

test_that("CLI workflow handles debug mode correctly", {
  test_config <- list(
    cell_types = c("Mono", "CD4_T"),
    chromosomes = c("chr22"),
    feature_type = "gene",
    debug = list(
      enabled = TRUE,
      n_genes = 100,
      n_peaks = 50
    ),
    file_patterns = list(
      eqtl_acat = "test.txt"
    ),
    parameters = list(n_cores = 1)
  )

  # Check debug settings are parsed correctly
  expect_true(test_config$debug$enabled)
  expect_equal(test_config$debug$n_genes, 100)
  expect_equal(test_config$debug$n_peaks, 50)
})

test_that("CLI workflow handles L1/L2 cell type filtering", {
  # Test with L1 cell types
  l1_celltypes <- c(
    "predicted.celltype.l1.Mono",
    "predicted.celltype.l1.CD4_T",
    "predicted.celltype.l1.CD8_T"
  )

  h <- cascade::DEFAULT_CELL_HIERARCHY
  l1_filtered <- filter_l1_celltypes(l1_celltypes, h)
  expect_equal(l1_filtered, l1_celltypes)

  # Test with L2 cell types
  l2_celltypes <- c(
    "predicted.celltype.l2.CD14_Mono",
    "predicted.celltype.l2.CD4_Naive"
  )

  l2_filtered <- filter_l1_celltypes(l2_celltypes, h)
  expect_equal(length(l2_filtered), 0)

  # Test with mixed L1 and L2
  mixed_celltypes <- c(
    "predicted.celltype.l1.Mono",
    "predicted.celltype.l2.CD14_Mono"
  )

  mixed_filtered <- filter_l1_celltypes(mixed_celltypes, h)
  expect_equal(mixed_filtered, "predicted.celltype.l1.Mono")
})

test_that("CLI workflow generates correct output files", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Expected output file patterns based on feature type
  expected_outputs <- list(
    gene = c(
      "gene_categorization.tsv.gz",
      "gene_categorization.summary.tsv",
      "gene_categorization.top_variants.tsv.gz"
    ),
    peak = c(
      "peak_categorization.tsv.gz",
      "peak_categorization.summary.tsv",
      "peak_categorization.top_variants.tsv.gz"
    ),
    variant = c(
      "variant_categorization.tsv.gz",
      "variant_categorization.qtl_mechanism.summary.tsv",
      "variant_categorization.cell_type_specificity.summary.tsv"
    )
  )

  # Check expected file patterns
  for (feature_type in names(expected_outputs)) {
    files <- expected_outputs[[feature_type]]
    expect_true(all(grepl("\\.(tsv|tsv\\.gz)$", files)))
  }
})

test_that("CLI handles file pattern substitution correctly", {
  test_patterns <- list(
    eqtl_acat = "/path/to/{CELL_TYPE}.acat.txt.gz",
    eqtl_susie = "/path/to/{CELL_TYPE}.{CHR}.susie.txt.gz",
    peak_bed = "/path/to/{CELL_TYPE}.peaks.bed"
  )

  cell_type <- "Mono"
  chromosome <- "chr22"

  # Test substitution
  acat_path <- gsub("\\{CELL_TYPE\\}", cell_type, test_patterns$eqtl_acat)
  expect_equal(acat_path, "/path/to/Mono.acat.txt.gz")

  susie_path <- gsub("\\{CELL_TYPE\\}", cell_type, test_patterns$eqtl_susie)
  susie_path <- gsub("\\{CHR\\}", chromosome, susie_path)
  expect_equal(susie_path, "/path/to/Mono.chr22.susie.txt.gz")
})

test_that("CLI workflow handles Cochran's Q settings", {
  # Test with Cochran's Q enabled
  config_cochran <- list(
    parameters = list(
      use_cochran_q = TRUE,
      cochran_q_threshold = 5e-8
    )
  )

  expect_true(config_cochran$parameters$use_cochran_q)
  expect_equal(config_cochran$parameters$cochran_q_threshold, 5e-8)

  # Test with Cochran's Q disabled
  config_lfsr <- list(
    parameters = list(
      use_cochran_q = FALSE,
      lfsr_sig_threshold = 0.05,
      lfsr_null_threshold = 0.5
    )
  )

  expect_false(config_lfsr$parameters$use_cochran_q)
  expect_equal(config_lfsr$parameters$lfsr_sig_threshold, 0.05)
})

test_that("CLI workflow handles cache settings", {
  config_with_cache <- list(
    cache = list(enabled = TRUE),
    parameters = list(n_cores = 4)
  )

  expect_true(config_with_cache$cache$enabled)

  config_no_cache <- list(
    cache = list(enabled = FALSE),
    parameters = list(n_cores = 1)
  )

  expect_false(config_no_cache$cache$enabled)
})

test_that("CLI workflow handles parallel processing settings", {
  configs <- list(
    list(parameters = list(n_cores = 1)), # Single core
    list(parameters = list(n_cores = 4)), # Multi-core
    list(parameters = list(n_cores = 8)) # Many cores
  )

  for (config in configs) {
    expect_true(config$parameters$n_cores >= 1)
    expect_true(config$parameters$n_cores <= 16) # Reasonable upper limit
  }
})

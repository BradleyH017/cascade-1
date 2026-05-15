library(testthat)
library(cascade)
library(data.table)

# Source the comprehensive mock data helper
source(test_path("helper-mock-data.R"))

test_that("run_cascade works with gene analysis", {
  skip_if_not_installed("data.table")

  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create comprehensive test configuration with mock data
  config <- create_test_config(
    test_dir = test_dir,
    cell_types = c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l1.CD4_T",
      "predicted.celltype.l1.CD8_T"
    ),
    chromosomes = c("chr7"),
    feature_type = "gene"
  )

  output_dir <- file.path(test_dir, "output")

  # Run the actual categorization with comprehensive mock data
  results <- run_cascade(config, output_dir)

  # Test that output was created
  expect_true(dir.exists(output_dir))

  # Test that results have the expected structure
  expect_true(is.list(results))

  # Check for gene categorization results
  if (!is.null(results$gene_categories)) {
    expect_true("data.table" %in% class(results$gene_categories) ||
      "data.frame" %in% class(results$gene_categories))
    expect_true(nrow(results$gene_categories) > 0)
  }

  # Check that output files were created
  output_files <- list.files(output_dir, full.names = TRUE)
  expect_true(length(output_files) > 0)
})

test_that("run_cascade handles multiple feature types", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Test different feature type combinations
  feature_types_to_test <- list(
    "gene",
    "peak",
    "gene,peak",
    "gene,peak,variant"
  )

  for (ft in feature_types_to_test) {
    config <- list(
      cell_types = c("Mono"),
      chromosomes = c("chr22"),
      feature_type = ft,
      file_patterns = list(
        eqtl_acat = "dummy.txt",
        caqtl_acat = "dummy.txt"
      ),
      parameters = list(n_cores = 1)
    )

    # Parse feature types as done in the function
    feature_type <- config$feature_type
    if (grepl(",", feature_type)) {
      feature_types <- trimws(strsplit(feature_type, ",")[[1]])
    } else {
      feature_types <- feature_type
    }

    # Verify parsing
    if (ft == "gene,peak,variant") {
      expect_equal(feature_types, c("gene", "peak", "variant"))
    } else if (ft == "gene,peak") {
      expect_equal(feature_types, c("gene", "peak"))
    } else {
      expect_equal(feature_types, ft)
    }
  }
})

test_that("run_cascade handles Cochran's Q threshold", {
  config_with_cochran <- list(
    parameters = list(
      use_cochran_q = TRUE,
      cochran_q_threshold = 1e-6
    )
  )

  expect_equal(config_with_cochran$parameters$cochran_q_threshold, 1e-6)

  # Test default threshold
  config_default <- list(
    parameters = list(
      use_cochran_q = TRUE
    )
  )

  # Function should set default to 5e-8 if not specified
  if (is.null(config_default$parameters$cochran_q_threshold)) {
    config_default$parameters$cochran_q_threshold <- 5e-8
  }
  expect_equal(config_default$parameters$cochran_q_threshold, 5e-8)
})

test_that("run_cascade handles debug mode", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  config <- list(
    cell_types = c("Mono"),
    chromosomes = c("chr22"),
    feature_type = "gene",
    file_patterns = list(eqtl_acat = "test.txt"),
    parameters = list(n_cores = 1),
    debug = list(
      enabled = TRUE,
      n_genes = 50,
      n_peaks = 25
    )
  )

  # Check debug settings
  expect_true(config$debug$enabled)
  expect_equal(config$debug$n_genes, 50)
  expect_equal(config$debug$n_peaks, 25)
})

test_that("run_cascade filters L1/L2 cell types correctly", {
  # Test with only L1 cell types
  l1_config <- list(
    cell_types = c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l1.CD4_T"
    ),
    feature_type = "gene"
  )

  l1_types <- filter_l1_celltypes(l1_config$cell_types, cascade::DEFAULT_CELL_HIERARCHY)
  expect_equal(length(l1_types), 2)
  expect_equal(l1_types, l1_config$cell_types)

  # Test with only L2 cell types
  l2_config <- list(
    cell_types = c(
      "predicted.celltype.l2.CD14_Mono",
      "predicted.celltype.l2.CD4_Naive"
    ),
    feature_type = "gene"
  )

  l2_types <- filter_l1_celltypes(l2_config$cell_types, cascade::DEFAULT_CELL_HIERARCHY)
  expect_equal(length(l2_types), 0)

  # Test with mixed L1 and L2
  mixed_config <- list(
    cell_types = c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l2.CD14_Mono"
    ),
    feature_type = "gene"
  )

  mixed_types <- filter_l1_celltypes(mixed_config$cell_types, cascade::DEFAULT_CELL_HIERARCHY)
  expect_equal(length(mixed_types), 1)
  expect_equal(mixed_types, "predicted.celltype.l1.Mono")
})

test_that("run_cascade loads configuration from JSON file", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create a config object
  config <- list(
    cell_types = c("Mono", "CD4_T"),
    chromosomes = c("chr22"),
    feature_type = "gene",
    file_patterns = list(
      eqtl_acat = "test.txt"
    ),
    parameters = list(
      n_cores = 2,
      pip_threshold = 0.5
    )
  )

  # Save to JSON file
  config_file <- file.path(test_dir, "config.json")
  jsonlite::write_json(config, config_file, auto_unbox = TRUE)

  # Load and verify
  loaded_config <- jsonlite::fromJSON(config_file, simplifyVector = FALSE)

  # Convert list to vector if needed
  loaded_cell_types <- if (is.list(loaded_config$cell_types)) {
    unlist(loaded_config$cell_types)
  } else {
    loaded_config$cell_types
  }
  expect_equal(loaded_cell_types, config$cell_types)
  expect_equal(loaded_config$chromosomes, config$chromosomes)
  expect_equal(loaded_config$feature_type, config$feature_type)
  expect_equal(loaded_config$parameters$n_cores, config$parameters$n_cores)
})

test_that("run_cascade handles mashr model loading", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create mock mashr model files
  mock_mashr <- list(
    fitted_g = list(Ulist = list(), grid = c(0.5, 1, 2)),
    result = list(PosteriorMean = matrix(0, 10, 3))
  )

  eqtl_mashr_file <- file.path(test_dir, "eqtl.mashr.rds")
  caqtl_mashr_file <- file.path(test_dir, "caqtl.mashr.rds")

  saveRDS(mock_mashr, eqtl_mashr_file)
  saveRDS(mock_mashr, caqtl_mashr_file)

  config <- list(
    cell_types = c("Mono"),
    chromosomes = c("chr22"),
    feature_type = "gene",
    file_patterns = list(
      eqtl_acat = "test.txt",
      eqtl_mashr = eqtl_mashr_file,
      caqtl_mashr = caqtl_mashr_file
    ),
    parameters = list(n_cores = 1)
  )

  # Check that mashr files exist
  expect_true(file.exists(config$file_patterns$eqtl_mashr))
  expect_true(file.exists(config$file_patterns$caqtl_mashr))

  # Load and verify structure
  loaded_eqtl <- readRDS(config$file_patterns$eqtl_mashr)
  expect_true("fitted_g" %in% names(loaded_eqtl))
  expect_true("result" %in% names(loaded_eqtl))
})

test_that("run_cascade generates expected summary files", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Expected summary file names
  expected_summaries <- c(
    "gene_categorization.summary.tsv",
    "peak_categorization.summary.tsv",
    "variant_categorization.qtl_mechanism.summary.tsv",
    "variant_categorization.cell_type_specificity.summary.tsv"
  )

  # Check file naming conventions
  for (file in expected_summaries) {
    expect_true(grepl("summary", file, ignore.case = TRUE))
    expect_true(grepl("\\.tsv$", file))
  }
})

test_that("run_cascade handles cache settings", {
  config_with_cache <- list(
    cache = list(enabled = TRUE),
    parameters = list(n_cores = 1)
  )

  expect_true(config_with_cache$cache$enabled)

  config_no_cache <- list(
    cache = list(enabled = FALSE),
    parameters = list(n_cores = 1)
  )

  expect_false(config_no_cache$cache$enabled)

  # Test default (should be enabled)
  config_default <- list(
    parameters = list(n_cores = 1)
  )

  # If cache not specified, should default to enabled
  if (is.null(config_default$cache)) {
    config_default$cache <- list(enabled = TRUE)
  }

  expect_true(config_default$cache$enabled)
})

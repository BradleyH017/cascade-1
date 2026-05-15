library(testthat)
library(cascade)
library(jsonlite)
library(data.table)

# Source the mock data helper
source(test_path("helper-mock-data.R"))

# These end-to-end tests invoke the installed CLI script via Rscript in a
# subprocess. Under R CMD check the package lives in a check-local library
# that the subprocess cannot always discover, so the subprocess fails
# silently (system(intern = TRUE) does not surface non-zero exit codes).
# Skip these on CRAN-style checks and CI matrix runs unless explicitly
# enabled via CASCADE_RUN_CLI_TESTS.
skip_if_not_cli_runnable <- function() {
  testthat::skip_on_cran()
  if (!nzchar(Sys.getenv("CASCADE_RUN_CLI_TESTS"))) {
    testthat::skip("Set CASCADE_RUN_CLI_TESTS=1 to enable CLI subprocess tests")
  }
}

test_that("CLI workflow runs end-to-end with gene analysis", {
  skip_if_not_cli_runnable()
  skip_if_not_installed("argparse")
  skip_if_not_installed("cli")

  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create test configuration with mock data
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

  # Save configuration to JSON file
  config_file <- file.path(test_dir, "test_config.json")
  write_json(config, config_file, auto_unbox = TRUE, pretty = TRUE)

  # Define output directory
  output_dir <- file.path(test_dir, "results")

  # Run the CLI script programmatically
  cli_script <- system.file("scripts", "cascade_cli.R", package = "cascade")

  if (cli_script != "" && file.exists(cli_script)) {
    # Build command
    cmd <- sprintf(
      "Rscript '%s' --config '%s' --output_dir '%s'",
      cli_script, config_file, output_dir
    )

    # Execute CLI command
    result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

    # Check that the command succeeded
    expect_true(dir.exists(output_dir))

    # Check for expected output files
    output_files <- list.files(output_dir,
      pattern = "\\.tsv|\\.txt",
      full.names = FALSE
    )
    expect_true(length(output_files) > 0)

    # Check for gene categorization output
    gene_cat_file <- file.path(output_dir, "gene_categorization.tsv.gz")
    if (file.exists(gene_cat_file)) {
      gene_data <- fread(gene_cat_file)
      expect_true(nrow(gene_data) > 0)
      # Check for gene identifier column (could be gene_id, feature, or gene)
      expect_true("gene_id" %in% names(gene_data) || "feature" %in% names(gene_data) || "gene" %in% names(gene_data))
    }
  } else {
    # If CLI script not installed, test direct function call
    results <- run_cascade(config, output_dir)
    expect_true(dir.exists(output_dir))
    expect_true(is.list(results))
  }
})

test_that("CLI workflow handles peak analysis", {
  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create test configuration for peaks
  config <- create_test_config(
    test_dir = test_dir,
    cell_types = c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l1.CD4_T"
    ),
    chromosomes = c("chr7"),
    feature_type = "peak"
  )

  output_dir <- file.path(test_dir, "results")

  # Run categorization
  results <- run_cascade(config, output_dir)

  # Check results
  expect_true(dir.exists(output_dir))

  # Check for peak-specific outputs
  if (!is.null(results$peak_categories)) {
    expect_true(nrow(results$peak_categories) > 0)
  }
})

test_that("CLI workflow handles variant analysis", {
  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create test configuration for variants
  config <- create_test_config(
    test_dir = test_dir,
    cell_types = c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l1.CD4_T",
      "predicted.celltype.l1.CD8_T"
    ),
    chromosomes = c("chr7"),
    feature_type = "variant"
  )

  output_dir <- file.path(test_dir, "results")

  # Run categorization
  results <- run_cascade(config, output_dir)

  # Check results
  expect_true(dir.exists(output_dir))

  # Check for variant-specific outputs
  if (!is.null(results$variant_categories)) {
    # variant_categories is a list with per_celltype and cross_celltype components
    expect_true(is.list(results$variant_categories))
    if (!is.null(results$variant_categories$cross_celltype)) {
      expect_true(nrow(results$variant_categories$cross_celltype) > 0)
    }
  }
})

test_that("CLI workflow handles all feature types", {
  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create test configuration for all features
  config <- create_test_config(
    test_dir = test_dir,
    cell_types = c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l1.CD4_T"
    ),
    chromosomes = c("chr7"),
    feature_type = "all"
  )

  output_dir <- file.path(test_dir, "results")

  # Run categorization
  results <- run_cascade(config, output_dir)

  # Check results
  expect_true(dir.exists(output_dir))
  expect_true(is.list(results))

  # Check that multiple feature types were processed
  features_processed <- c(
    !is.null(results$gene_categories),
    !is.null(results$peak_categories),
    !is.null(results$variant_categories)
  )
  expect_true(sum(features_processed) >= 2)
})

test_that("CLI workflow respects debug mode", {
  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create configuration with debug mode
  config <- create_test_config(
    test_dir = test_dir,
    cell_types = c("predicted.celltype.l1.Mono"),
    chromosomes = c("chr7"),
    feature_type = "gene"
  )

  # Enable debug mode with small limits
  config$debug <- list(
    enabled = TRUE,
    n_genes = 10,
    n_peaks = 5
  )

  output_dir <- file.path(test_dir, "results")

  # Run with debug mode
  results <- run_cascade(config, output_dir)

  # Check that debug mode limited processing
  if (!is.null(results$gene_categories)) {
    # In debug mode, should process limited number of genes
    expect_true(nrow(results$gene_categories) <= config$debug$n_genes * 2)
  }
})

test_that("CLI workflow handles L1 and L2 cell types correctly", {
  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Test with mixed L1 and L2 cell types
  all_cell_types <- c(
    "predicted.celltype.l1.Mono",
    "predicted.celltype.l1.CD4_T",
    "predicted.celltype.l2.CD14_Mono",
    "predicted.celltype.l2.CD4_Naive"
  )

  config <- create_test_config(
    test_dir = test_dir,
    cell_types = all_cell_types,
    chromosomes = c("chr7"),
    feature_type = "gene"
  )

  output_dir <- file.path(test_dir, "results")

  # Run categorization
  results <- run_cascade(config, output_dir)

  # The function should filter to use only L1 cell types
  expect_true(dir.exists(output_dir))

  # Check that L1 filtering was applied
  l1_types <- filter_l1_celltypes(all_cell_types, cascade::DEFAULT_CELL_HIERARCHY)
  expect_equal(length(l1_types), 2) # Should only have 2 L1 types
})

test_that("CLI workflow handles Cochran's Q settings", {
  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Test with Cochran's Q enabled
  config <- create_test_config(
    test_dir = test_dir,
    cell_types = c("predicted.celltype.l1.Mono"),
    chromosomes = c("chr7"),
    feature_type = "gene"
  )

  # Explicitly set Cochran's Q parameters
  config$parameters$use_cochran_q <- TRUE
  config$parameters$cochran_q_threshold <- 1e-6

  output_dir <- file.path(test_dir, "results")

  # Run categorization
  results <- run_cascade(config, output_dir)

  # Check that results were produced
  expect_true(dir.exists(output_dir))
})

test_that("CLI workflow saves and loads JSON configuration correctly", {
  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create configuration
  config <- create_test_config(
    test_dir = test_dir,
    cell_types = c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l1.CD4_T"
    ),
    chromosomes = c("chr7", "chr22"),
    feature_type = "gene"
  )

  # Save to JSON
  config_file <- file.path(test_dir, "config.json")
  write_json(config, config_file, auto_unbox = TRUE, pretty = TRUE)

  # Reload and verify
  loaded_config <- fromJSON(config_file, simplifyVector = FALSE)

  # Convert cell types to vector if needed
  loaded_cell_types <- if (is.list(loaded_config$cell_types)) {
    unlist(loaded_config$cell_types)
  } else {
    loaded_config$cell_types
  }

  expect_equal(loaded_cell_types, config$cell_types)

  # Convert chromosomes to vector if needed
  loaded_chromosomes <- if (is.list(loaded_config$chromosomes)) {
    unlist(loaded_config$chromosomes)
  } else {
    loaded_config$chromosomes
  }
  expect_equal(loaded_chromosomes, config$chromosomes)
  expect_equal(loaded_config$feature_type, config$feature_type)
  expect_equal(
    loaded_config$parameters$pip_threshold,
    config$parameters$pip_threshold
  )
})

test_that("CLI workflow generates summary files", {
  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create configuration
  config <- create_test_config(
    test_dir = test_dir,
    cell_types = c("predicted.celltype.l1.Mono"),
    chromosomes = c("chr7"),
    feature_type = "gene"
  )

  output_dir <- file.path(test_dir, "results")

  # Run categorization
  results <- run_cascade(config, output_dir)

  # Check for summary files
  summary_files <- list.files(output_dir,
    pattern = "summary",
    full.names = FALSE
  )

  # Should have at least one summary file
  expect_true(length(summary_files) >= 1)
})

test_that("CLI workflow handles multiple chromosomes", {
  test_dir <- tempfile("cascade_cli_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create configuration with multiple chromosomes
  config <- create_test_config(
    test_dir = test_dir,
    cell_types = c("predicted.celltype.l1.Mono"),
    chromosomes = c("chr7", "chr22"),
    feature_type = "gene"
  )

  output_dir <- file.path(test_dir, "results")

  # Run categorization
  results <- run_cascade(config, output_dir)

  # Check results
  expect_true(dir.exists(output_dir))
  expect_true(is.list(results))
})

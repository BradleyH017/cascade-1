#!/usr/bin/env Rscript

#' Command Line Interface for Cascade Package
#' 
#' Provides a command-line interface for running cascade analysis

suppressPackageStartupMessages({
  if (!requireNamespace("argparse", quietly = TRUE)) {
    stop("Package 'argparse' is required. Install it with: install.packages('argparse')")
  }
  library(argparse)
  library(cascade)
  library(jsonlite)
  library(cli)
})

# Create argument parser
parser <- ArgumentParser(
  description = "CASCADE: Cell type specificity categorization for QTL analysis",
  epilog = "For more information, see the CASCADE documentation"
)

# Add arguments
parser$add_argument(
  "--config",
  type = "character",
  required = TRUE,
  help = "Configuration file in JSON format"
)

parser$add_argument(
  "--output_dir",
  type = "character",
  default = "results",
  help = "Output directory for results [default: %(default)s]"
)

parser$add_argument(
  "--feature_type",
  type = "character",
  choices = c("gene", "peak", "variant", "all"),
  help = "Feature type to analyze (overrides config file)"
)

parser$add_argument(
  "--version",
  action = "store_true",
  help = "Show package version and exit"
)

# Parse arguments
args <- parser$parse_args()

# Show version if requested
if (args$version) {
  cli::cli_alert_info("CASCADE version: {packageVersion('cascade')}")
  quit(status = 0)
}

# Validate config file
if (!file.exists(args$config)) {
  cli::cli_abort("Configuration file not found: {.file {args$config}}")
}

# Load configuration
cli::cli_h1("CASCADE Analysis")
cli::cli_alert_info("Loading configuration from {.file {args$config}}")

config <- tryCatch(
  jsonlite::fromJSON(args$config, simplifyVector = FALSE),
  error = function(e) {
    cli::cli_abort("Failed to parse configuration file: {e$message}")
  }
)

# Override feature type if provided
if (!is.null(args$feature_type)) {
  config$feature_type <- args$feature_type
  cli::cli_alert_info("Overriding feature type to: {.val {args$feature_type}}")
}

# Display analysis parameters
cli::cli_h2("Analysis Parameters")
cli::cli_ul(c(
  "Output directory: {.path {args$output_dir}}",
  "Feature type: {.val {config$feature_type}}",
  "Cell types: {.val {paste(config$cell_types, collapse = ', ')}}",
  "Chromosomes: {.val {paste(config$chromosomes, collapse = ', ')}}"
))

# Check debug mode
if (!is.null(config$debug) && config$debug$enabled) {
  cli::cli_alert_warning("DEBUG MODE ENABLED")
  cli::cli_ul(c(
    "Processing up to {.val {config$debug$n_genes}} genes",
    "Processing up to {.val {config$debug$n_peaks}} peaks"
  ))
}

# Run analysis
cli::cli_h2("Running Analysis")
cli::cli_progress_step("Initializing CASCADE analysis...")

tryCatch({
  # Create output directory
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Run the analysis
  results <- run_cascade(config, args$output_dir)
  
  # Display results summary
  cli::cli_h2("Results Summary")
  
  if (!is.null(results$gene_categories)) {
    cli::cli_alert_success("Gene categorization: {.val {nrow(results$gene_categories)}} genes processed")
  }
  
  if (!is.null(results$peak_categories)) {
    cli::cli_alert_success("Peak categorization: {.val {nrow(results$peak_categories)}} peaks processed")
  }
  
  if (!is.null(results$variant_categories)) {
    cli::cli_alert_success("Variant categorization: {.val {nrow(results$variant_categories)}} variants processed")
  }
  
  # List output files
  output_files <- list.files(args$output_dir, pattern = "\\.(tsv|txt|rds|tsv\\.gz)$", full.names = TRUE)
  if (length(output_files) > 0) {
    cli::cli_h3("Output Files")
    for (f in output_files) {
      size <- file.info(f)$size
      size_str <- if (size > 1e6) {
        sprintf("%.1f MB", size / 1e6)
      } else if (size > 1e3) {
        sprintf("%.1f KB", size / 1e3)
      } else {
        sprintf("%d bytes", size)
      }
      cli::cli_alert_info("{.file {basename(f)}} ({size_str})")
    }
  }
  
  cli::cli_alert_success("Analysis completed successfully!")
  
}, error = function(e) {
  cli::cli_alert_danger("Error during analysis: {e$message}")
  
  # Provide helpful error messages for common issues
  if (grepl("No ACAT files", e$message)) {
    cli::cli_alert_info("Check that your file patterns in the config match your data files")
  } else if (grepl("subscript out of bounds", e$message)) {
    cli::cli_alert_info("This may indicate empty or malformed input files")
  }
  
  quit(status = 1)
})
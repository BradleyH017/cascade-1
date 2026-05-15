#' Two-stage variant categorization
#'
#' Main entry point for variant categorization. Runs per-cell-type QTL pattern
#' detection (Stage 1), then aggregates across cell types with bulk
#' `data.table` operations (Stage 2).
#'
#' @param variant_data Variant data from `load_variant_data_by_qtl_type()`.
#' @param lfsr_results Pre-loaded LFSR results (optional).
#' @param output_dir Optional output directory for saving intermediate results.
#' @param config Configuration object (optional).
#' @return List with `per_celltype`, `cross_celltype`, and
#'   `variant_feature_specificity` results.
#' @export
categorize_variants <- function(variant_data, lfsr_results = NULL, output_dir = NULL, config = NULL) {
  n_variants <- length(variant_data$all_variants)
  cli::cli_alert_info("Categorizing {.val {n_variants}} variants")

  # Set mc.cores for parallel Stage 1 processing if config provides n_cores
  if (!is.null(config$parameters$n_cores) && config$parameters$n_cores > 1) {
    old_mc_cores <- getOption("mc.cores", 1L)
    options(mc.cores = as.integer(config$parameters$n_cores))
    on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  }

  cli::cli_h3("Stage 1: Per-Cell-Type Variant Categorization")

  # Prepare the variant matrix
  cli::cli_alert_info("Preparing variant matrix...")
  variant_matrix <- prepare_variant_matrix(
    variant_data = variant_data,
    lfsr_results = lfsr_results
  )
  cli::cli_alert_success("Variant matrix prepared with {nrow(variant_matrix)} rows")

  cli::cli_alert_info("Preparing matrix structures...")
  matrices <- prepare_variant_matrices(variant_matrix, variant_data)

  n_matrix_variants <- length(matrices$variant_ids)
  cli::cli_alert_info("Processing {n_matrix_variants} variants across {length(matrices$cell_types)} cell types")

  # Detect if we have L2 cell types
  all_cell_types <- matrices$cell_types
  hierarchy <- if (!is.null(config$hierarchy)) config$hierarchy else DEFAULT_CELL_HIERARCHY
  l2_cell_types <- all_cell_types[is_l2_celltype(all_cell_types, hierarchy)]
  l1_cell_types <- filter_l1_celltypes(all_cell_types, hierarchy)

  if (length(l2_cell_types) > 0) {
    cli::cli_alert_info("Detected {length(l2_cell_types)} L2 cell type(s) - will aggregate together with L1 cell types")
    cli::cli_alert_info("Cross-cell-type aggregation will use all {length(all_cell_types)} cell type(s) (L1 + L2)")
  }

  # Get all peak-gene links (will be filtered by cell type later)
  all_peak_gene_links <- if (!is.null(variant_data$peak_gene_links) && nrow(variant_data$peak_gene_links) > 0) {
    data.table::as.data.table(variant_data$peak_gene_links)
  } else {
    data.table::data.table(peak_id = character(0), gene_id = character(0))
  }

  # Stage 1: Process each cell type separately (both L1 and L2)
  per_celltype_results <- process_celltype_variants(
    cell_types = matrices$cell_types,
    matrices = matrices,
    all_peak_gene_links = all_peak_gene_links,
    output_dir = output_dir
  )

  cli::cli_h3("Stage 2: Cross-Cell-Type Variant Categorization (Combined L1 + L2)")

  # Combined L1 + L2 aggregation
  if (length(all_cell_types) > 0) {
    cli::cli_alert_info("Aggregating all cell types ({length(all_cell_types)} types: {length(l1_cell_types)} L1 + {length(l2_cell_types)} L2)")

    results_with_specificity <- aggregate_variant_categorization(
      variant_ids = matrices$variant_ids,
      cell_types = all_cell_types,
      per_celltype_results = per_celltype_results,
      lfsr_results = lfsr_results,
      config = config
    )

    combined <- results_with_specificity$categorization
    variant_feature_specificity <- results_with_specificity$variant_feature_specificity
  } else {
    cli::cli_alert_warning("No cell types available for cross-cell-type aggregation")
    combined <- data.frame(
      variant_id = character(0),
      qtl_mechanism_category = character(0),
      qtl_pattern_number = integer(0),
      affected_cell_types = character(0),
      n_affected_cell_types = integer(0),
      best_cell_type = character(0),
      significant_cell_types = character(0),
      qtl_patterns = character(0),
      cell_type_specificity = character(0),
      associated_genes = character(0),
      associated_peaks = character(0),
      cascade_peak_genes = character(0),
      peak_overlap = logical(0),
      caqtl_status = character(0),
      peak_gene_link = character(0),
      eqtl_status = character(0),
      stringsAsFactors = FALSE
    )
    variant_feature_specificity <- NULL
  }

  cli::cli_alert_success("Bulk categorization completed successfully")

  # Save cross-cell-type results if output directory specified
  if (!is.null(output_dir)) {
    if (!is.null(combined) && nrow(combined) > 0) {
      combined_file <- file.path(output_dir, "variant_categorization.tsv")
      utils::write.table(combined, combined_file,
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      cli::cli_alert_success("Saved combined L1+L2 cross-cell-type results to: {.file {combined_file}}")
    }
  }

  return(list(
    per_celltype = per_celltype_results,
    cross_celltype = combined,
    variant_feature_specificity = variant_feature_specificity
  ))
}

#' Save Per-Cell-Type Variant Results
#'
#' Save Stage 1 results to separate TSV files per cell type
#'
#' @param per_celltype_results List of per-cell-type results
#' @param output_dir Output directory
#' @import data.table
save_variant_results_per_celltype <- function(per_celltype_results, output_dir) {
  for (ct in names(per_celltype_results)) {
    ct_results <- per_celltype_results[[ct]]

    # Create filename following design specification
    output_file <- file.path(
      output_dir,
      paste0("variant_categorization_per_celltype_", ct, ".tsv.gz")
    )

    # Write results
    data.table::fwrite(ct_results, output_file,
      sep = "\t", quote = FALSE, compress = "gzip", na = "NA"
    )

    cli::cli_alert_success("Saved {ct} results to: {.file {output_file}}")
  }
}

#' Save Cross-Cell-Type Variant Results
#'
#' Save Stage 2 results to TSV file (combined L1 + L2)
#'
#' @param cross_celltype_results Cross-cell-type results data frame
#' @param output_dir Output directory
#' @param suffix Optional suffix for file names (e.g., "l2")
#' @import data.table
save_variant_results_cross_celltype <- function(cross_celltype_results, output_dir, suffix = "") {
  prefix <- if (suffix != "") paste0("variant_categorization_", suffix) else "variant_categorization"
  output_file <- file.path(output_dir, paste0(prefix, ".tsv.gz"))
  qtl_summary_file <- file.path(output_dir, paste0(prefix, ".qtl_mechanism.summary.tsv"))
  spec_summary_file <- file.path(output_dir, paste0(prefix, ".cell_type_specificity.summary.tsv"))

  data.table::fwrite(cross_celltype_results, output_file,
    sep = "\t", quote = FALSE, compress = "gzip", na = "NA"
  )

  cli::cli_alert_success("Saved cross-cell-type results to: {.file {output_file}}")

  qtl_summary <- summarize_column(cross_celltype_results, "qtl_mechanism_category", "QTL Mechanism")
  save_categorization_summary(qtl_summary, qtl_summary_file)

  specificity_summary <- summarize_column(cross_celltype_results, "cell_type_specificity", "Cell Type Specificity")
  save_categorization_summary(specificity_summary, spec_summary_file)
}

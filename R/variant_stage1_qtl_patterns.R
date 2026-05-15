#' @title Table Categorization Stage 1
#' @description Per-cell-type QTL pattern detection functions
#' Process Cell Type Variants
#'
#' Performs Stage 1 per-cell-type variant categorization
#' @param cell_types Vector of cell type names
#' @param matrices Prepared variant matrices
#' @param all_peak_gene_links Peak-gene links data.table
#' @param output_dir Optional output directory
#' @return List of per-cell-type results
#' @noRd
process_celltype_variants <- function(cell_types, matrices, all_peak_gene_links, output_dir = NULL) {
  # Create variant matrix data for C++ function
  variant_matrix_data <- list(
    variant_ids = matrices$variant_ids,
    cell_types = matrices$cell_types,
    eqtl_matrix = matrices$variant_matrix_wide$n_eqtl_matrix > 0,
    caqtl_matrix = matrices$variant_matrix_wide$n_caqtl_matrix > 0,
    caqtl_peaks_list = matrices$caqtl_peaks_list,
    eqtl_genes_list = matrices$eqtl_genes_list,
    peak_overlaps = matrices$peak_overlaps_ordered
  )

  # Pre-compute constant arguments
  pattern_to_mechanism <- sapply(QTL_PATTERNS, function(p) p$mechanism)

  # Process each cell type - use parallel if available and beneficial
  n_cores <- getOption("mc.cores", 1L)
  use_parallel <- .Platform$OS.type == "unix" && n_cores > 1 && length(cell_types) > 1

  process_one_ct <- function(ct_idx) {
    ct <- cell_types[ct_idx]

    # Filter peak-gene links for this cell type
    ct_peak_gene_links <- filter_peak_gene_links_for_celltype(
      all_peak_gene_links, ct
    )

    # Call the existing C++ pattern detection for this cell type
    ct_results <- detect_variant_qtl_patterns(
      variant_ids = matrices$variant_ids,
      variant_matrix_data = variant_matrix_data,
      peak_gene_links = ct_peak_gene_links,
      cell_types = matrices$cell_types,
      ct_idx = as.integer(ct_idx - 1), # 0-based index for C++
      qtl_mechanism_categories = QTL_MECHANISMS,
      pattern_to_mechanism = pattern_to_mechanism,
      caqtl_status_values = CAQTL_STATUS,
      eqtl_status_values = EQTL_STATUS
    )

    ct_results
  }

  if (use_parallel) {
    # Parallel processing: mclapply uses fork(), so shared data is copy-on-write
    cli::cli_alert_info("Using parallel processing for Stage 1 ({n_cores} cores)")
    result_list <- parallel::mclapply(
      seq_along(cell_types),
      process_one_ct,
      mc.cores = min(n_cores, length(cell_types))
    )
    # Check for errors from mclapply
    errors <- vapply(result_list, inherits, logical(1), "try-error")
    if (any(errors)) {
      cli::cli_alert_warning("Parallel processing failed for {sum(errors)} cell type(s), falling back to sequential")
      for (i in which(errors)) {
        result_list[[i]] <- process_one_ct(i)
      }
    }
  } else {
    # Sequential processing
    result_list <- vector("list", length(cell_types))
    for (ct_idx in seq_along(cell_types)) {
      cli::cli_alert("Processing {cell_types[ct_idx]}")
      result_list[[ct_idx]] <- process_one_ct(ct_idx)
    }
  }

  # Name results by cell type
  per_celltype_results <- setNames(result_list, cell_types)

  # Optionally save per-cell-type results
  if (!is.null(output_dir)) {
    for (ct in cell_types) {
      ct_file <- file.path(output_dir, paste0("variant_categorization_per_celltype_", ct, ".tsv"))
      write.table(per_celltype_results[[ct]], ct_file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  return(per_celltype_results)
}

#' Filter Peak-Gene Links for Cell Type
#'
#' Filters peak-gene links for a specific cell type
#' @param all_peak_gene_links Data.table with all peak-gene links
#' @param target_cell_type Cell type to filter for
#' @return Filtered data.table
#' @noRd
filter_peak_gene_links_for_celltype <- function(all_peak_gene_links, target_cell_type) {
  # Select columns: always peak_id and gene_id, plus mechanism if available
  has_mechanism <- "mechanism" %in% names(all_peak_gene_links)
  select_cols <- if (has_mechanism) c("peak_id", "gene_id", "mechanism") else c("peak_id", "gene_id")

  if (nrow(all_peak_gene_links) > 0) {
    if (!"config_cell_type" %in% names(all_peak_gene_links)) {
      cli::cli_abort("Peak-gene links missing required {.field config_cell_type} column. Available: {.val {names(all_peak_gene_links)}}")
    }
    ct_peak_gene_links <- all_peak_gene_links[config_cell_type == target_cell_type, ..select_cols]
  } else {
    ct_peak_gene_links <- data.table::data.table(peak_id = character(0), gene_id = character(0))
    if (has_mechanism) ct_peak_gene_links[, mechanism := character(0)]
  }
  return(ct_peak_gene_links)
}

#' Categorize Variants Using Table Format
#'
#' Performs variant categorization using a clean variant × cell_type table
#' with proper two-stage analysis: per-cell-type pattern detection followed
#' by cross-cell-type aggregation
#'
#' @param variant_data Variant data from load_variant_data() (must include peak_gene_links)
#' @param lfsr_results LFSR results from load_lfsr_results()
#' @param output_dir Optional directory to save per-cell-type results
#' @return DataFrame with cross-cell-type categorization results
#' @export

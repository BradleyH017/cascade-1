#' @title Table Categorization Utilities
#' @description Matrix preparation and helper functions for variant categorization
#' Prepare Variant Matrices for Categorization
#'
#' Converts variant data into matrix form for vectorized C++ processing.
#' @param variant_matrix Prepared variant matrix from prepare_variant_matrix()
#' @param variant_data Original variant data including peak overlaps
#' @return List containing variant matrices and metadata
#' @noRd
prepare_variant_matrices <- function(variant_matrix, variant_data) {
  # Get unique variants and cell types
  variant_ids <- unique(variant_matrix$variant_id)
  cell_types <- unique(variant_matrix$cell_type)
  n_variants <- length(variant_ids)

  # Convert to data.table for efficient operations
  dt <- data.table::as.data.table(variant_matrix)

  # Create index mappings for fast lookup
  variant_idx_map <- setNames(seq_along(variant_ids), variant_ids)
  celltype_idx_map <- setNames(seq_along(cell_types), cell_types)

  # Add indices to data.table
  dt[, `:=`(
    var_idx = variant_idx_map[variant_id],
    ct_idx = celltype_idx_map[cell_type]
  )]

  # Initialize matrices with proper dimensions
  n_eqtl_matrix <- matrix(0, nrow = n_variants, ncol = length(cell_types))
  n_caqtl_matrix <- matrix(0, nrow = n_variants, ncol = length(cell_types))
  eqtl_genes_matrix <- matrix("", nrow = n_variants, ncol = length(cell_types))
  caqtl_peaks_matrix <- matrix("", nrow = n_variants, ncol = length(cell_types))

  # Vectorized matrix filling using data.table indexing
  matrix_indices <- as.matrix(dt[, .(var_idx, ct_idx)])

  n_eqtl_matrix[matrix_indices] <- dt$n_eqtl_genes
  n_caqtl_matrix[matrix_indices] <- dt$n_caqtl_peaks

  # Handle NA values for gene/peak strings
  dt[, `:=`(
    eqtl_genes = ifelse(is.na(eqtl_genes), "", eqtl_genes),
    caqtl_peaks = ifelse(is.na(caqtl_peaks), "", caqtl_peaks)
  )]

  eqtl_genes_matrix[matrix_indices] <- dt$eqtl_genes
  caqtl_peaks_matrix[matrix_indices] <- dt$caqtl_peaks

  # Create variant matrix wide format
  variant_matrix_wide <- list(
    variant_ids = variant_ids,
    cell_types = cell_types,
    n_eqtl_matrix = n_eqtl_matrix,
    n_caqtl_matrix = n_caqtl_matrix,
    eqtl_genes_matrix = eqtl_genes_matrix,
    caqtl_peaks_matrix = caqtl_peaks_matrix
  )

  # Create QTL lists from matrices
  qtl_lists <- create_qtl_lists_from_matrices(
    dt, n_variants, cell_types, variant_idx_map
  )

  # Process peak overlaps
  peak_overlaps_ordered <- process_peak_overlaps(
    variant_data$variant_peak_overlaps,
    variant_ids
  )

  # Combine all data
  return(list(
    variant_ids = variant_ids,
    cell_types = cell_types,
    variant_idx_map = variant_idx_map,
    variant_matrix_wide = variant_matrix_wide,
    caqtl_peaks_list = qtl_lists$caqtl_peaks_list,
    eqtl_genes_list = qtl_lists$eqtl_genes_list,
    peak_overlaps_ordered = peak_overlaps_ordered
  ))
}

#' Create QTL Lists from Variant Matrix
#'
#' Efficiently creates structured lists of QTL peaks and genes
#' @param dt Data.table with variant matrix data
#' @param n_variants Number of variants
#' @param cell_types Vector of cell type names
#' @param variant_idx_map Mapping of variant IDs to indices
#' @return List with caqtl_peaks_list and eqtl_genes_list
#' @noRd
create_qtl_lists_from_matrices <- function(dt, n_variants, cell_types, variant_idx_map) {
  # Pre-allocate the full list structure
  total_elements <- n_variants * length(cell_types)
  caqtl_peaks_list <- vector("list", total_elements)
  eqtl_genes_list <- vector("list", total_elements)

  # Create a mapping of matrix position to list position
  dt[, list_idx := var_idx + (ct_idx - 1) * n_variants]

  # Process non-empty peaks
  peaks_data <- dt[caqtl_peaks != "", .(list_idx, caqtl_peaks)]
  if (nrow(peaks_data) > 0) {
    peaks_split <- strsplit(peaks_data$caqtl_peaks, ",", fixed = TRUE)
    caqtl_peaks_list[peaks_data$list_idx] <- peaks_split
  }

  # Process non-empty genes
  genes_data <- dt[eqtl_genes != "", .(list_idx, eqtl_genes)]
  if (nrow(genes_data) > 0) {
    genes_split <- strsplit(genes_data$eqtl_genes, ",", fixed = TRUE)
    eqtl_genes_list[genes_data$list_idx] <- genes_split
  }

  # Fill remaining NULL positions with empty character vectors
  null_indices <- which(sapply(caqtl_peaks_list, is.null))
  caqtl_peaks_list[null_indices] <- list(character(0))

  null_indices <- which(sapply(eqtl_genes_list, is.null))
  eqtl_genes_list[null_indices] <- list(character(0))

  # Set dimensions for compatibility
  dim(caqtl_peaks_list) <- c(n_variants, length(cell_types))
  dim(eqtl_genes_list) <- c(n_variants, length(cell_types))

  return(list(
    caqtl_peaks_list = caqtl_peaks_list,
    eqtl_genes_list = eqtl_genes_list
  ))
}

#' Process Peak Overlaps for Variants
#'
#' Orders peak overlaps according to variant ID ordering
#' @param variant_peak_overlaps Peak overlaps from variant data
#' @param variant_ids Ordered variant IDs
#' @return Ordered list of peak overlaps
#' @noRd
process_peak_overlaps <- function(variant_peak_overlaps, variant_ids) {
  if (!is.null(variant_peak_overlaps) && length(variant_peak_overlaps) > 0) {
    # Use match for vectorized lookup - directly index by variant IDs
    peak_overlaps_ordered <- variant_peak_overlaps[variant_ids]
    # Replace NA/NULL with NULL explicitly for consistency
    peak_overlaps_ordered[sapply(peak_overlaps_ordered, is.null)] <- list(NULL)
  } else {
    peak_overlaps_ordered <- vector("list", length(variant_ids))
  }
  return(peak_overlaps_ordered)
}

#' Build Pattern and Category Matrices
#'
#' Creates matrices of QTL patterns and categories from per-cell-type results
#' @param variant_ids Vector of variant IDs
#' @param cell_types Vector of cell type names
#' @param per_celltype_results List of per-cell-type results
#' @return List with pattern_matrix, category_matrix, variant_idx_map
#' @noRd
build_pattern_category_matrices <- function(variant_ids, cell_types, per_celltype_results) {
  n_variants <- length(variant_ids)

  # Initialize matrices
  pattern_matrix <- matrix(25, nrow = n_variants, ncol = length(cell_types))
  category_matrix <- matrix(QTL_MECHANISMS[8], nrow = n_variants, ncol = length(cell_types))
  colnames(pattern_matrix) <- cell_types
  colnames(category_matrix) <- cell_types

  # Create variant index lookup for O(1) access
  variant_idx_map <- setNames(seq_along(variant_ids), variant_ids)

  for (ct_idx in seq_along(cell_types)) {
    ct <- cell_types[ct_idx]
    if (ct %in% names(per_celltype_results)) {
      ct_results <- per_celltype_results[[ct]]

      # Use vectorized matching
      matching_indices <- variant_idx_map[ct_results$variant_id]
      valid_matches <- !is.na(matching_indices)

      if (any(valid_matches)) {
        pattern_matrix[matching_indices[valid_matches], ct_idx] <- ct_results$qtl_pattern_number[valid_matches]
        category_matrix[matching_indices[valid_matches], ct_idx] <- ct_results$qtl_mechanism_category[valid_matches]
      }
    }
  }

  return(list(
    pattern_matrix = pattern_matrix,
    category_matrix = category_matrix,
    variant_idx_map = variant_idx_map
  ))
}

#' Aggregate Variant Features from Cell Type Results
#'
#' Aggregates genes, peaks, and cascade information across cell types
#' @param n_variants Number of variants
#' @param per_celltype_results List of per-cell-type results
#' @param variant_idx_map Mapping of variant IDs to indices
#' @return List with aggregated variant features
#' @noRd
aggregate_variant_features <- function(n_variants, per_celltype_results, variant_idx_map) {
  # Combine all cell type results into one data.table
  all_results <- data.table::rbindlist(
    lapply(names(per_celltype_results), function(ct) {
      dt <- data.table::as.data.table(per_celltype_results[[ct]])
      dt[, cell_type := ct]
      dt
    }),
    fill = TRUE
  )

  # Add variant indices for efficient mapping
  all_results[, var_idx := variant_idx_map[variant_id]]

  # Process genes using data.table aggregation
  gene_agg <- process_feature_column(
    all_results,
    "associated_genes",
    variant_idx_map,
    n_variants
  )

  # Process peaks using data.table aggregation
  peak_agg <- process_feature_column(
    all_results,
    "associated_peaks",
    variant_idx_map,
    n_variants
  )

  # Process cascade peak-gene pairs
  cascade_agg <- process_cascade_pairs(
    all_results,
    variant_idx_map,
    n_variants
  )

  return(list(
    variant_genes_map = gene_agg$feature_map,
    variant_peaks_map = peak_agg$feature_map,
    variant_associated_genes = gene_agg$associated_vector,
    variant_associated_peaks = peak_agg$associated_vector,
    variant_cascade_peak_genes = cascade_agg$cascade_vector
  ))
}

#' Process Feature Column using data.table
#'
#' Vectorized helper to process gene or peak columns during aggregation
#' @noRd
process_feature_column <- function(all_results, column_name, variant_idx_map, n_variants) {
  # Filter to valid features
  dt <- all_results[
    !is.na(get(column_name)) & get(column_name) != "" & get(column_name) != "NA",
    .(variant_id, var_idx, features = get(column_name))
  ]

  if (nrow(dt) == 0) {
    return(list(
      feature_map = list(),
      associated_vector = character(n_variants)
    ))
  }

  # Split features and unnest
  dt[, features_list := strsplit(features, ",", fixed = TRUE)]
  dt_unnested <- dt[, .(feature = unlist(features_list)), by = .(variant_id, var_idx)]

  # Create feature map by aggregating unique features per variant
  feature_map_dt <- dt_unnested[, .(features = list(unique(feature))), by = variant_id]
  feature_map <- as.list(setNames(feature_map_dt$features, feature_map_dt$variant_id))

  # Create associated vector (comma-separated unique features per variant)
  assoc_dt <- dt_unnested[, .(associated = paste(unique(feature), collapse = ",")),
    by = var_idx
  ]

  # Initialize full vector with empty strings
  associated_vector <- character(n_variants)
  # Fill in the values at the appropriate indices
  associated_vector[assoc_dt$var_idx] <- assoc_dt$associated

  return(list(
    feature_map = feature_map,
    associated_vector = associated_vector
  ))
}

#' Process Cascade Peak-Gene Pairs using data.table
#'
#' Vectorized helper to process cascade peak-gene pairs during aggregation
#' @noRd
process_cascade_pairs <- function(all_results, variant_idx_map, n_variants) {
  # Filter to valid cascade pairs
  dt <- all_results[
    !is.na(cascade_peak_genes) & cascade_peak_genes != "" & cascade_peak_genes != "NA",
    .(variant_id, var_idx, pairs = cascade_peak_genes)
  ]

  if (nrow(dt) == 0) {
    return(list(
      cascade_vector = character(n_variants)
    ))
  }

  # Split pairs and unnest
  dt[, pairs_list := strsplit(pairs, ",", fixed = TRUE)]
  dt_unnested <- dt[, .(pair = unlist(pairs_list)), by = .(variant_id, var_idx)]

  # Aggregate unique pairs per variant
  cascade_dt <- dt_unnested[, .(cascade = paste(unique(pair), collapse = ",")),
    by = var_idx
  ]

  # Initialize full vector with empty strings
  cascade_vector <- character(n_variants)
  # Fill in the values at the appropriate indices
  cascade_vector[cascade_dt$var_idx] <- cascade_dt$cascade

  return(list(
    cascade_vector = cascade_vector
  ))
}

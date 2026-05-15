#' @title Table Categorization Stage 2
#' @description Cross-cell-type aggregation using bulk data.table operations
#' Process QTL Variant Results
#'
#' Processes results for variants with QTL effects
#' @param qtl_variants Indices of QTL variants
#' @param variant_ids Vector of variant IDs
#' @param pattern_matrix Pattern matrix
#' @param category_matrix Category matrix
#' @param has_qtl_matrix Boolean matrix of QTL presence
#' @param cell_types Vector of cell types
#' @param per_celltype_results Per-cell-type results
#' @return List with processed QTL results
#' @noRd
process_qtl_variant_results <- function(qtl_variants, variant_ids, pattern_matrix,
                                        category_matrix, has_qtl_matrix, has_eqtl_matrix,
                                        has_caqtl_matrix, cell_types,
                                        per_celltype_results, best_patterns) {
  n_variants <- length(variant_ids)
  n_ct <- length(cell_types)

  # Initialize result vectors with defaults
  result_mechanisms <- character(n_variants)
  result_pattern_numbers <- integer(n_variants)
  result_affected_cts <- character(n_variants)
  result_gene_affected_cts <- character(n_variants)
  result_peak_affected_cts <- character(n_variants)
  variant_best_cts <- character(n_variants)
  variant_significant_cts <- character(n_variants)
  variant_qtl_patterns <- rep("No Effect", n_variants)
  variant_peak_overlap <- logical(n_variants)
  variant_caqtl_status <- rep("no_caqtl", n_variants)
  variant_peak_gene_link <- character(n_variants)
  variant_eqtl_status <- rep("no_eqtl", n_variants)
  variant_link_mechanism <- rep(NA_character_, n_variants)

  # --- Vectorized: pattern numbers and best CT index ---
  result_pattern_numbers[qtl_variants] <- best_patterns[qtl_variants]

  # For each QTL variant, find the first column matching best_pattern
  # best_ct_idx_vec[i] = column index of first cell type with best pattern for variant i
  best_ct_idx_vec <- integer(n_variants)
  for (i in qtl_variants) {
    best_ct_idx_vec[i] <- which(pattern_matrix[i, ] == best_patterns[i])[1]
  }

  # Vectorized mechanism lookup using matrix indexing
  qtl_matrix_idx <- cbind(qtl_variants, best_ct_idx_vec[qtl_variants])
  result_mechanisms[qtl_variants] <- category_matrix[qtl_matrix_idx]

  # --- Vectorized: comma-separated cell type strings ---
  # Use data.table long-format grouping instead of apply() on 119K × 33 matrices.
  # Helper: convert boolean matrix to comma-separated CT strings per variant
  bool_matrix_to_ct_strings <- function(bool_mat, var_indices) {
    # Get (row, col) pairs where TRUE
    rc <- which(bool_mat[var_indices, , drop = FALSE], arr.ind = TRUE)
    if (nrow(rc) == 0) {
      return(rep("", length(var_indices)))
    }
    dt <- data.table::data.table(
      orig_row = rc[, 1], # row index within var_indices subset
      ct = cell_types[rc[, 2]]
    )
    agg <- dt[, .(cts = paste(ct, collapse = ",")), keyby = orig_row]
    result <- rep("", length(var_indices))
    result[agg$orig_row] <- agg$cts
    result
  }

  variant_significant_cts[qtl_variants] <- bool_matrix_to_ct_strings(has_qtl_matrix, qtl_variants)
  result_gene_affected_cts[qtl_variants] <- bool_matrix_to_ct_strings(has_eqtl_matrix, qtl_variants)
  result_peak_affected_cts[qtl_variants] <- bool_matrix_to_ct_strings(has_caqtl_matrix, qtl_variants)

  # Best CTs: pattern matches best_pattern (vectorized comparison)
  best_pattern_match <- pattern_matrix == best_patterns # recycled column-wise
  variant_best_cts[qtl_variants] <- bool_matrix_to_ct_strings(best_pattern_match, qtl_variants)

  # Affected CTs: best mechanism + has QTL
  best_mechs <- category_matrix[cbind(seq_len(n_variants), best_ct_idx_vec)]
  best_mech_match <- category_matrix == best_mechs & has_qtl_matrix # recycled
  result_affected_cts[qtl_variants] <- bool_matrix_to_ct_strings(best_mech_match, qtl_variants)

  # --- Bulk lookup: per-celltype result fields for best CT ---
  # Flatten all per-celltype results into one data.table for a single keyed join
  required_ct_cols <- c(
    "variant_id", "qtl_pattern_number", "peak_overlap",
    "caqtl", "peak_gene_link", "eqtl", "link_mechanism"
  )
  ct_data_list <- lapply(cell_types, function(ct) {
    ct_data <- per_celltype_results[[ct]]
    if (is.null(ct_data) || length(ct_data$variant_id) == 0) {
      return(NULL)
    }
    missing <- setdiff(required_ct_cols, names(ct_data))
    if (length(missing) > 0) {
      cli::cli_abort("Per-celltype results for {ct} missing required columns: {.val {missing}}")
    }
    vids <- ct_data$variant_id
    if (length(vids) == 0) {
      return(NULL)
    }
    data.table::data.table(
      variant_id = vids,
      cell_type = ct,
      qtl_pattern_number = ct_data$qtl_pattern_number,
      peak_overlap = ct_data$peak_overlap,
      caqtl = ct_data$caqtl,
      peak_gene_link = ct_data$peak_gene_link,
      eqtl = ct_data$eqtl,
      link_mechanism = ct_data$link_mechanism
    )
  })
  all_ct_dt <- data.table::rbindlist(ct_data_list[!vapply(ct_data_list, is.null, logical(1))])
  data.table::setkey(all_ct_dt, variant_id, cell_type)

  # Build lookup table: for each QTL variant, the best cell type
  best_ct_names <- cell_types[best_ct_idx_vec[qtl_variants]]
  lookup_dt <- data.table::data.table(
    variant_id = variant_ids[qtl_variants],
    cell_type = best_ct_names,
    orig_idx = qtl_variants
  )

  # Single keyed join for best-CT fields
  matched <- all_ct_dt[lookup_dt, on = .(variant_id, cell_type), nomatch = NA]

  # Fill results from matched rows (use caqtl as existence check — always present in C++ output)
  valid <- !is.na(matched$caqtl)
  if (any(valid)) {
    idx <- matched$orig_idx[valid]
    variant_peak_overlap[idx] <- matched$peak_overlap[valid]
    variant_caqtl_status[idx] <- matched$caqtl[valid]
    variant_peak_gene_link[idx] <- matched$peak_gene_link[valid]
    variant_eqtl_status[idx] <- matched$eqtl[valid]

    # Map qtl_pattern_number to human-readable pattern interpretation
    # (Fixes bug: old code looked for non-existent "qtl_pattern" field,
    #  so qtl_pattern was always "No Effect" even for QTL variants)
    pnums <- matched$qtl_pattern_number[valid]
    pattern_interp <- vapply(QTL_PATTERNS, function(p) p$interpretation, character(1))
    valid_pnum <- !is.na(pnums) & pnums >= 1 & pnums <= length(pattern_interp)
    variant_qtl_patterns[idx[valid_pnum]] <- pattern_interp[pnums[valid_pnum]]
  }

  # --- Bulk link_mechanism aggregation ---
  # For each QTL variant, aggregate link_mechanism across all significant CTs
  # Build a long-format table of (variant_id, significant_cell_type) pairs
  sig_pairs_list <- lapply(qtl_variants, function(i) {
    sig_ct <- cell_types[has_qtl_matrix[i, ]]
    if (length(sig_ct) > 0) {
      data.table::data.table(variant_id = variant_ids[i], cell_type = sig_ct, orig_idx = i)
    } else {
      NULL
    }
  })
  sig_pairs <- data.table::rbindlist(sig_pairs_list[!vapply(sig_pairs_list, is.null, logical(1))])

  if (nrow(sig_pairs) > 0) {
    # Join with flattened per-CT data to get link_mechanism
    sig_mechs <- all_ct_dt[sig_pairs,
      on = .(variant_id, cell_type),
      .(variant_id, orig_idx, link_mechanism)
    ]
    # Filter to non-empty mechanisms
    sig_mechs <- sig_mechs[!is.na(link_mechanism) & nchar(link_mechanism) > 0]

    if (nrow(sig_mechs) > 0) {
      # Split comma-separated mechanisms, aggregate unique sorted per variant
      mech_agg <- sig_mechs[,
        {
          all_mechs <- unlist(strsplit(link_mechanism, ","))
          .(link_mechanism = paste(sort(unique(all_mechs)), collapse = ","))
        },
        by = .(orig_idx)
      ]
      variant_link_mechanism[mech_agg$orig_idx] <- mech_agg$link_mechanism
    }
  }

  return(list(
    result_mechanisms = result_mechanisms,
    result_pattern_numbers = result_pattern_numbers,
    result_affected_cts = result_affected_cts,
    result_gene_affected_cts = result_gene_affected_cts,
    result_peak_affected_cts = result_peak_affected_cts,
    variant_best_cts = variant_best_cts,
    variant_significant_cts = variant_significant_cts,
    variant_qtl_patterns = variant_qtl_patterns,
    variant_peak_overlap = variant_peak_overlap,
    variant_caqtl_status = variant_caqtl_status,
    variant_peak_gene_link = variant_peak_gene_link,
    variant_eqtl_status = variant_eqtl_status,
    variant_link_mechanism = variant_link_mechanism
  ))
}

#' Stage 2 Cross-Cell-Type Variant Categorization
#'
#' Uses bulk data.table operations for LFSR lookups.
#'
#' @param variant_ids Vector of variant IDs
#' @param cell_types Vector of cell type names
#' @param per_celltype_results List of per-cell-type results
#' @param lfsr_results LFSR results object
#' @return Data frame with cross-cell-type categorization results
#' @noRd
aggregate_variant_categorization <- function(variant_ids, cell_types,
                                             per_celltype_results,
                                             lfsr_results = NULL,
                                             config = NULL) {
  n_variants <- length(variant_ids)

  # Extract hierarchy-derived C++ params (or use defaults)
  hierarchy <- if (!is.null(config$hierarchy)) config$hierarchy else DEFAULT_CELL_HIERARCHY
  cpp_params <- hierarchy_to_cpp_params(hierarchy)

  if (n_variants > 100) {
    cli::cli_alert_info("Processing {n_variants} variants using bulk data.table approach")
  }

  # Step 1: Build pattern and category matrices
  matrices <- build_pattern_category_matrices(variant_ids, cell_types, per_celltype_results)
  pattern_matrix <- matrices$pattern_matrix
  category_matrix <- matrices$category_matrix
  variant_idx_map <- matrices$variant_idx_map

  # Step 2: Vectorized pattern analysis
  # Use direct matrix %in% instead of element-wise apply() to avoid per-cell function call overhead.
  no_qtl_patterns <- c(23L, 24L, 25L)
  has_qtl_matrix <- matrix(!(pattern_matrix %in% no_qtl_patterns),
    nrow = nrow(pattern_matrix), ncol = ncol(pattern_matrix)
  )

  # Find best patterns and significant cell types
  best_patterns <- do.call(pmin, as.data.frame(pattern_matrix))
  n_affected_cts <- rowSums(has_qtl_matrix)

  # Step 2b: Separate eQTL and caQTL affected cell types
  eqtl_patterns <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 20L, 21L, 22L)
  caqtl_patterns <- c(1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L)

  has_eqtl_matrix <- matrix(pattern_matrix %in% eqtl_patterns,
    nrow = nrow(pattern_matrix), ncol = ncol(pattern_matrix)
  )
  has_caqtl_matrix <- matrix(pattern_matrix %in% caqtl_patterns,
    nrow = nrow(pattern_matrix), ncol = ncol(pattern_matrix)
  )

  n_gene_affected_cts <- rowSums(has_eqtl_matrix)
  n_peak_affected_cts <- rowSums(has_caqtl_matrix)

  # Step 3: Pre-extract all variant-gene-peak mappings in bulk
  cli::cli_alert_info("Pre-extracting variant-feature mappings (bulk operation)...")

  features <- aggregate_variant_features(n_variants, per_celltype_results, variant_idx_map)
  variant_genes_map <- features$variant_genes_map
  variant_peaks_map <- features$variant_peaks_map
  variant_associated_genes <- features$variant_associated_genes
  variant_associated_peaks <- features$variant_associated_peaks
  variant_cascade_peak_genes <- features$variant_cascade_peak_genes

  # Step 4: Bulk LFSR processing using data.table merge operations
  cli::cli_alert_info("Processing LFSR values using bulk data.table operations...")

  # Prepare variant features for bulk lookup
  variant_features <- list(
    variant_ids = variant_ids,
    variant_genes_map = variant_genes_map,
    variant_peaks_map = variant_peaks_map
  )

  # Perform bulk LFSR lookup - THIS IS THE KEY OPTIMIZATION
  if (!is.null(lfsr_results)) {
    # Single bulk operation replaces all individual lookups
    lfsr_lookup_results <- bulk_variant_lfsr_lookup(
      variant_features = variant_features,
      lfsr_results = lfsr_results,
      cell_types = cell_types
    )
  } else {
    # Create dummy LFSR results if no data
    lfsr_lookup_results <- data.table::data.table(
      variant_id = rep(variant_ids, each = length(cell_types)),
      cell_type = rep(cell_types, length(variant_ids)),
      min_lfsr = 1.0
    )
  }

  # Step 5: Process results
  cli::cli_alert_info("Finalizing categorization results...")

  # Process in vectorized batches
  no_qtl_variants <- which(n_affected_cts == 0)
  qtl_variants <- which(n_affected_cts > 0)

  # Initialize result vectors
  result_mechanisms <- character(n_variants)
  result_pattern_numbers <- integer(n_variants)
  result_affected_cts <- character(n_variants)

  # Handle no-QTL variants in batch
  if (length(no_qtl_variants) > 0) {
    result_mechanisms[no_qtl_variants] <- QTL_MECHANISMS[8]
    result_pattern_numbers[no_qtl_variants] <- 25
    result_affected_cts[no_qtl_variants] <- ""
  }

  # Process QTL variants
  qtl_results <- if (length(qtl_variants) > 0) {
    process_qtl_variant_results(
      qtl_variants, variant_ids, pattern_matrix, category_matrix,
      has_qtl_matrix, has_eqtl_matrix, has_caqtl_matrix, cell_types,
      per_celltype_results, best_patterns
    )
  } else {
    list(
      result_mechanisms = character(n_variants),
      result_pattern_numbers = integer(n_variants),
      result_affected_cts = character(n_variants),
      result_gene_affected_cts = character(n_variants),
      result_peak_affected_cts = character(n_variants),
      variant_best_cts = character(n_variants),
      variant_significant_cts = character(n_variants),
      variant_qtl_patterns = character(n_variants),
      variant_peak_overlap = logical(n_variants),
      variant_caqtl_status = character(n_variants),
      variant_peak_gene_link = character(n_variants),
      variant_eqtl_status = character(n_variants),
      variant_link_mechanism = rep(NA_character_, n_variants)
    )
  }

  # Merge QTL results
  if (length(qtl_variants) > 0) {
    result_mechanisms[qtl_variants] <- qtl_results$result_mechanisms[qtl_variants]
    result_pattern_numbers[qtl_variants] <- qtl_results$result_pattern_numbers[qtl_variants]
    result_affected_cts[qtl_variants] <- qtl_results$result_affected_cts[qtl_variants]
  }

  # Create preliminary results for bulk specificity categorization
  preliminary_results <- data.table::data.table(
    variant_id = variant_ids,
    qtl_mechanism_category = result_mechanisms,
    n_affected_cell_types = n_affected_cts,
    n_gene_affected_cell_types = n_gene_affected_cts,
    n_peak_affected_cell_types = n_peak_affected_cts,
    affected_cell_types = result_affected_cts,
    gene_affected_cell_types = if (length(qtl_variants) > 0) qtl_results$result_gene_affected_cts else character(n_variants),
    peak_affected_cell_types = if (length(qtl_variants) > 0) qtl_results$result_peak_affected_cts else character(n_variants),
    qtl_pattern_number = result_pattern_numbers,
    qtl_pattern = {
      # Build qtl_pattern from pattern numbers: QTL variants get interpretations from
      # process_qtl_variant_results; no-QTL variants get mapped here from their pattern number
      qp <- if (length(qtl_variants) > 0) qtl_results$variant_qtl_patterns else character(n_variants)
      if (length(no_qtl_variants) > 0) {
        no_qtl_interp <- vapply(
          result_pattern_numbers[no_qtl_variants],
          get_pattern_interpretation, character(1)
        )
        qp[no_qtl_variants] <- no_qtl_interp
      }
      qp
    },
    best_cell_types = if (length(qtl_variants) > 0) qtl_results$variant_best_cts else character(n_variants),
    significant_cts = if (length(qtl_variants) > 0) qtl_results$variant_significant_cts else character(n_variants),
    associated_genes = variant_associated_genes,
    associated_peaks = variant_associated_peaks,
    cascade_peak_genes = variant_cascade_peak_genes,
    peak_overlap = if (length(qtl_variants) > 0) qtl_results$variant_peak_overlap else logical(n_variants),
    caqtl = if (length(qtl_variants) > 0) qtl_results$variant_caqtl_status else character(n_variants),
    peak_gene_link = if (length(qtl_variants) > 0) qtl_results$variant_peak_gene_link else character(n_variants),
    eqtl = if (length(qtl_variants) > 0) qtl_results$variant_eqtl_status else character(n_variants),
    link_mechanism = if (length(qtl_variants) > 0) qtl_results$variant_link_mechanism else rep(NA_character_, n_variants)
  )

  # Bulk categorize cell type specificity
  cli::cli_alert_info("Categorizing cell type specificity (bulk operation)...")

  # Combined specificity using mechanism-based approach without LFSR checking
  # Pass NULL for LFSR to disable "Likely shared but underpowered" detection
  specificities <- categorize_variants_bulk_cpp(
    variant_results = preliminary_results,
    lfsr_lookup_results = NULL, # No LFSR checking for combined mechanism-based specificity
    l2_to_l1_mapping = cpp_params$l2_to_l1_mapping,
    lineage_groups = cpp_params$lineage_groups,
    subgroup_levels = cpp_params$subgroup_levels,
    bulk_cts = cpp_params$bulk_cts,
    other_cts = cpp_params$other_cts,
    lfsr_sig_threshold = config$parameters$lfsr_sig_threshold,
    lfsr_null_threshold = config$parameters$lfsr_null_threshold,
    specificity_categories = cpp_params$specificity_categories
  )

  # Gene-specific specificity: prepare data with gene-specific columns
  gene_results <- data.table::copy(preliminary_results)
  gene_results[, affected_cell_types := gene_affected_cell_types]
  gene_results[, n_affected_cell_types := n_gene_affected_cell_types]

  # Filter LFSR to only gene-specific values
  gene_lfsr <- if (!is.null(lfsr_lookup_results)) {
    lfsr_lookup_results[!is.na(min_lfsr_gene), .(variant_id, cell_type, min_lfsr = min_lfsr_gene)]
  } else {
    NULL
  }

  gene_specificities <- categorize_variants_bulk_cpp(
    variant_results = gene_results,
    lfsr_lookup_results = gene_lfsr,
    l2_to_l1_mapping = cpp_params$l2_to_l1_mapping,
    lineage_groups = cpp_params$lineage_groups,
    subgroup_levels = cpp_params$subgroup_levels,
    bulk_cts = cpp_params$bulk_cts,
    other_cts = cpp_params$other_cts,
    lfsr_sig_threshold = config$parameters$lfsr_sig_threshold,
    lfsr_null_threshold = config$parameters$lfsr_null_threshold,
    specificity_categories = cpp_params$specificity_categories
  )

  # Peak-specific specificity: prepare data with peak-specific columns
  peak_results <- data.table::copy(preliminary_results)
  peak_results[, affected_cell_types := peak_affected_cell_types]
  peak_results[, n_affected_cell_types := n_peak_affected_cell_types]

  # Filter LFSR to only peak-specific values
  peak_lfsr <- if (!is.null(lfsr_lookup_results)) {
    lfsr_lookup_results[!is.na(min_lfsr_peak), .(variant_id, cell_type, min_lfsr = min_lfsr_peak)]
  } else {
    NULL
  }

  peak_specificities <- categorize_variants_bulk_cpp(
    variant_results = peak_results,
    lfsr_lookup_results = peak_lfsr,
    l2_to_l1_mapping = cpp_params$l2_to_l1_mapping,
    lineage_groups = cpp_params$lineage_groups,
    subgroup_levels = cpp_params$subgroup_levels,
    bulk_cts = cpp_params$bulk_cts,
    other_cts = cpp_params$other_cts,
    lfsr_sig_threshold = config$parameters$lfsr_sig_threshold,
    lfsr_null_threshold = config$parameters$lfsr_null_threshold,
    specificity_categories = cpp_params$specificity_categories
  )

  # Add specificities to results
  preliminary_results$cell_type_specificity <- specificities
  preliminary_results$gene_cell_type_specificity <- gene_specificities
  preliminary_results$peak_cell_type_specificity <- peak_specificities

  # Reorder columns to match expected output format (exclude internal columns)
  final_results <- preliminary_results[, .(
    variant_id,
    qtl_mechanism_category,
    cell_type_specificity,
    gene_cell_type_specificity,
    peak_cell_type_specificity,
    qtl_pattern_number,
    qtl_pattern,
    best_cell_types,
    significant_cts,
    gene_affected_cell_types,
    peak_affected_cell_types,
    associated_genes,
    associated_peaks,
    cascade_peak_genes,
    peak_overlap,
    caqtl,
    peak_gene_link,
    eqtl,
    link_mechanism
  )]

  cli::cli_alert_success("Bulk Stage 2 completed for {n_variants} variants")

  # Compute variant-feature-level specificity using data.table joins + bulk C++ categorization.
  cli::cli_alert_info("Computing variant-feature-level specificity...")

  # Step 1: Build variant-feature long table from maps
  gene_pairs <- create_variant_feature_pairs(variant_ids, variant_genes_map, "gene")
  peak_pairs <- create_variant_feature_pairs(variant_ids, variant_peaks_map, "peak")

  vf_parts <- list()

  # Step 2: For each feature type, join with preliminary_results and categorize
  for (ft_info in list(
    list(
      pairs = gene_pairs, affected_col = "gene_affected_cell_types",
      n_affected_col = "n_gene_affected_cell_types",
      lfsr_col = "min_lfsr_gene", ft = "gene"
    ),
    list(
      pairs = peak_pairs, affected_col = "peak_affected_cell_types",
      n_affected_col = "n_peak_affected_cell_types",
      lfsr_col = "min_lfsr_peak", ft = "peak"
    )
  )) {
    if (nrow(ft_info$pairs) == 0) next

    # Join pairs with preliminary_results to get affected_cell_types
    pr_cols <- preliminary_results[, .(variant_id)]
    pr_cols[, affected_cell_types := preliminary_results[[ft_info$affected_col]]]
    pr_cols[, n_affected_cell_types := preliminary_results[[ft_info$n_affected_col]]]
    vf_dt <- merge(ft_info$pairs[, .(variant_id, feature_id)], pr_cols,
      by = "variant_id", all.x = FALSE
    )
    # Filter to variants that have affected cell types
    vf_dt <- vf_dt[!is.na(affected_cell_types) & nchar(affected_cell_types) > 0]

    if (nrow(vf_dt) == 0) next

    # Get feature-type-specific LFSR
    ft_lfsr <- if (!is.null(lfsr_lookup_results) && ft_info$lfsr_col %in% names(lfsr_lookup_results)) {
      lfsr_lookup_results[
        !is.na(get(ft_info$lfsr_col)),
        .(variant_id, cell_type, min_lfsr = get(ft_info$lfsr_col))
      ]
    } else {
      NULL
    }

    # Categorize using existing bulk C++ (one call for all pairs of this type)
    specificities <- categorize_variants_bulk_cpp(
      variant_results = vf_dt,
      lfsr_lookup_results = ft_lfsr,
      l2_to_l1_mapping = cpp_params$l2_to_l1_mapping,
      lineage_groups = cpp_params$lineage_groups,
      subgroup_levels = cpp_params$subgroup_levels,
      bulk_cts = cpp_params$bulk_cts,
      other_cts = cpp_params$other_cts,
      lfsr_sig_threshold = config$parameters$lfsr_sig_threshold,
      lfsr_null_threshold = config$parameters$lfsr_null_threshold,
      specificity_categories = cpp_params$specificity_categories
    )

    vf_dt[, `:=`(
      feature_type = ft_info$ft,
      variant_specificity = specificities,
      significant_cell_types = affected_cell_types
    )]

    vf_parts[[ft_info$ft]] <- vf_dt[, .(
      variant_id, feature_id, feature_type,
      variant_specificity, significant_cell_types
    )]
  }

  variant_feature_specificity <- data.table::rbindlist(vf_parts)
  cli::cli_alert_success("Computed variant-feature specificity for {nrow(variant_feature_specificity)} variant-feature pairs")

  return(list(
    categorization = final_results,
    variant_feature_specificity = variant_feature_specificity
  ))
}

#' Internal bulk LFSR lookup function
#' @noRd
bulk_variant_lfsr_lookup <- function(variant_features, lfsr_results, cell_types) {
  if (is.null(lfsr_results)) {
    # Return default high LFSR values if no data available
    return(data.table::data.table(
      variant_id = rep(variant_features$variant_ids, each = length(cell_types)),
      cell_type = rep(cell_types, length(variant_features$variant_ids)),
      min_lfsr = 1.0
    ))
  }

  # Step 1: Create variant-feature pairs for all variants
  variant_gene_pairs <- create_variant_feature_pairs(
    variant_features$variant_ids,
    variant_features$variant_genes_map,
    "gene"
  )

  variant_peak_pairs <- create_variant_feature_pairs(
    variant_features$variant_ids,
    variant_features$variant_peaks_map,
    "peak"
  )

  # Step 2: Use pre-computed long format from load_lfsr_results
  eqtl_lfsr_long <- if (!is.null(lfsr_results$eqtl_long) && nrow(variant_gene_pairs) > 0) {
    long_dt <- lfsr_results$eqtl_long
    if (length(cell_types) > 0) long_dt <- long_dt[cell_type %in% cell_types]
    long_dt
  } else {
    NULL
  }

  caqtl_lfsr_long <- if (!is.null(lfsr_results$caqtl_long) && nrow(variant_peak_pairs) > 0) {
    long_dt <- lfsr_results$caqtl_long
    if (length(cell_types) > 0) long_dt <- long_dt[cell_type %in% cell_types]
    long_dt
  } else {
    NULL
  }

  # Step 3: Perform bulk merges and aggregation (no CJ — build result from merges)
  gene_lfsr_results <- perform_lfsr_merge_and_aggregate(
    eqtl_lfsr_long,
    variant_gene_pairs,
    cell_types,
    result_col = "min_lfsr_gene"
  )

  peak_lfsr_results <- perform_lfsr_merge_and_aggregate(
    caqtl_lfsr_long,
    variant_peak_pairs,
    cell_types,
    result_col = "min_lfsr_peak"
  )

  # Step 4: Combine gene and peak LFSR values
  final_results <- combine_lfsr_results(
    variant_features$variant_ids,
    cell_types,
    gene_lfsr_results,
    peak_lfsr_results
  )

  return(final_results)
}

#' Helper: Create variant-feature pairs
#' @noRd
create_variant_feature_pairs <- function(variant_ids, variant_feature_map, feature_type) {
  # Vectorized: extract all features at once, build vectors, single DT construction
  # (replaces lapply + rbindlist of 119K tiny data.tables)
  non_empty <- vapply(variant_feature_map[variant_ids], length, integer(1)) > 0
  if (!any(non_empty)) {
    return(data.table::data.table(
      variant_id = character(0),
      feature_id = character(0),
      feature_type = character(0)
    ))
  }
  vids_with_features <- variant_ids[non_empty]
  feat_lists <- variant_feature_map[vids_with_features]
  feat_lens <- lengths(feat_lists)
  data.table::data.table(
    variant_id = rep(vids_with_features, feat_lens),
    feature_id = unlist(feat_lists, use.names = FALSE),
    feature_type = feature_type
  )
}

#' Helper: Perform LFSR merge and aggregation
#' @noRd
perform_lfsr_merge_and_aggregate <- function(lfsr_long, variant_feature_pairs,
                                             cell_types, result_col) {
  if (is.null(lfsr_long) || nrow(variant_feature_pairs) == 0) {
    return(NULL)
  }

  # Join variant-feature pairs directly with LFSR long data (no cross-join with cell types)
  # The LFSR long data already has (feature_id, variant_id, cell_type, lfsr)
  # We just need to find LFSR values for our variant-feature pairs
  data.table::setkey(variant_feature_pairs, feature_id, variant_id)

  # Semi-join: keep only LFSR rows that match our variant-feature pairs
  lfsr_matched <- lfsr_long[variant_feature_pairs,
    on = .(feature_id, variant_id),
    nomatch = NULL,
    .(variant_id, cell_type, lfsr)
  ]

  if (nrow(lfsr_matched) == 0) {
    return(NULL)
  }

  # Aggregate to get minimum LFSR per variant-cell type
  lfsr_results <- lfsr_matched[
    !is.na(lfsr),
    stats::setNames(list(min(lfsr)), result_col),
    by = .(variant_id, cell_type)
  ]
  data.table::setkey(lfsr_results, variant_id, cell_type)

  return(lfsr_results)
}

#' Helper: Combine gene and peak LFSR results
#' @noRd
combine_lfsr_results <- function(variant_ids, cell_types, gene_lfsr_results, peak_lfsr_results) {
  # Build result from LFSR data directly (avoid 119K × 33 CJ)
  # Only create rows that have at least one LFSR value
  if (!is.null(gene_lfsr_results) && !is.null(peak_lfsr_results)) {
    final_results <- merge(gene_lfsr_results, peak_lfsr_results,
      by = c("variant_id", "cell_type"), all = TRUE
    )
  } else if (!is.null(gene_lfsr_results)) {
    final_results <- data.table::copy(gene_lfsr_results)
    final_results[, min_lfsr_peak := NA_real_]
  } else if (!is.null(peak_lfsr_results)) {
    final_results <- data.table::copy(peak_lfsr_results)
    final_results[, min_lfsr_gene := NA_real_]
  } else {
    return(data.table::data.table(
      variant_id = character(0), cell_type = character(0),
      min_lfsr = NA_real_, min_lfsr_gene = NA_real_, min_lfsr_peak = NA_real_
    ))
  }

  # Take minimum across gene and peak LFSR
  if (!("min_lfsr_gene" %in% names(final_results))) final_results[, min_lfsr_gene := NA_real_]
  if (!("min_lfsr_peak" %in% names(final_results))) final_results[, min_lfsr_peak := NA_real_]
  final_results[, min_lfsr := pmin(min_lfsr_gene, min_lfsr_peak, na.rm = TRUE)]

  # Handle cases with no LFSR data
  final_results[is.na(min_lfsr), min_lfsr := 1.0]

  # Keep min_lfsr_gene and min_lfsr_peak for separate specificity categorization
  # They will be used for gene-specific and peak-specific specificity analysis

  return(final_results)
}

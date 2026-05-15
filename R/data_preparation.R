#' Prepare Feature Matrix for Categorization
#'
#' Creates a feature × cell_type matrix with ACAT q-values and minimum LFSR
#'
#' @param feature_data Loaded feature data with ACAT results
#' @param lfsr_results Loaded LFSR results
#' @param feature_type Either "gene" or "peak"
#' @return data.table with columns: feature_id, cell_type, acat_qval, min_lfsr
#' @importFrom data.table melt setDT setorder
prepare_feature_matrix <- function(feature_data, lfsr_results, feature_type = "gene") {
  # Extract ACAT matrix
  acat_df <- feature_data$acat_matrix
  feature_ids <- acat_df$feature_id
  cell_types <- setdiff(names(acat_df), "feature_id")

  # Convert to long format: feature_id × cell_type
  acat_long <- data.table::melt(data.table::setDT(acat_df),
    id.vars = "feature_id",
    measure.vars = cell_types,
    variable.name = "cell_type",
    value.name = "acat_qval"
  )

  # Get pre-computed feature-level LFSR lookup from new format
  lfsr_feature <- if (feature_type == "gene") {
    lfsr_results$eqtl_feature
  } else {
    lfsr_results$caqtl_feature
  }

  if (!is.null(lfsr_feature) && nrow(lfsr_feature) > 0) {
    # Merge with ACAT data using pre-computed feature-level LFSR
    data.table::setkeyv(acat_long, c("feature_id", "cell_type"))
    result <- merge(acat_long, lfsr_feature,
      by = c("feature_id", "cell_type"),
      all.x = TRUE
    )
  } else {
    result <- acat_long
    result[, min_lfsr := NA_real_]
  }

  # Ensure cell_type is character, not factor
  result[, cell_type := as.character(cell_type)]

  data.table::setorder(result, feature_id, cell_type)
  return(result)
}

#' Prepare Variant Matrix for Categorization
#'
#' Creates a variant by cell-type matrix with QTL effects and LFSR values
#'
#' @param variant_data Loaded variant data with SuSiE 95\% CS results (potentially filtered by PIP thresholds)
#' @param lfsr_results Loaded LFSR results
#' @return data.table with columns: variant_id, cell_type, eqtl_genes, caqtl_peaks, eqtl_lfsr, caqtl_lfsr
#' @noRd
prepare_variant_matrix <- function(variant_data, lfsr_results) {
  # Get unique variants and cell types
  if (is.null(variant_data$all_variants)) stop("No 'all_variants' field found in variant_data")
  variant_ids <- variant_data$all_variants

  cell_types <- if (!is.null(variant_data$data_by_ct)) {
    names(variant_data$data_by_ct)
  } else if (!is.null(variant_data$caqtl_data_by_ct)) {
    unique(c(names(variant_data$caqtl_data_by_ct), names(variant_data$eqtl_data_by_ct)))
  } else {
    stop("No cell type data found in variant_data")
  }

  # Create base matrix
  base_matrix <- data.table::CJ(variant_id = variant_ids, cell_type = cell_types)

  # Extract QTL effects from variant data
  qtl_effects <- data.table::rbindlist(lapply(cell_types, function(ct) {
    # Handle both old and new data structures
    ct_caqtl <- if (!is.null(variant_data$data_by_ct)) {
      variant_data$data_by_ct[[ct]]$caqtl
    } else if (!is.null(variant_data$caqtl_data_by_ct)) {
      variant_data$caqtl_data_by_ct[[ct]]
    } else {
      NULL
    }

    ct_eqtl <- if (!is.null(variant_data$data_by_ct)) {
      variant_data$data_by_ct[[ct]]$eqtl
    } else if (!is.null(variant_data$eqtl_data_by_ct)) {
      variant_data$eqtl_data_by_ct[[ct]]
    } else {
      NULL
    }

    # Get caQTL effects (data already has standardized column names from loader)
    caqtl_effects <- extract_qtl_effects(
      ct_caqtl,
      feature_col_out = "caqtl_peaks",
      count_col_out = "n_caqtl_peaks"
    )

    # Get eQTL effects (data already has standardized column names from loader)
    eqtl_effects <- extract_qtl_effects(
      ct_eqtl,
      feature_col_out = "eqtl_genes",
      count_col_out = "n_eqtl_genes"
    )

    # Combine for this cell type
    if (!is.null(caqtl_effects) || !is.null(eqtl_effects)) {
      # Start with variant IDs from either effect type
      variant_ids_ct <- unique(c(
        if (!is.null(caqtl_effects)) caqtl_effects$variant_id else character(0),
        if (!is.null(eqtl_effects)) eqtl_effects$variant_id else character(0)
      ))

      ct_result <- data.table::data.table(variant_id = variant_ids_ct, cell_type = ct)

      # Merge caQTL effects
      if (!is.null(caqtl_effects)) {
        ct_result <- merge(ct_result, caqtl_effects, by = "variant_id", all.x = TRUE)
      } else {
        ct_result[, ":="(caqtl_peaks = NA_character_, n_caqtl_peaks = 0L)]
      }

      # Merge eQTL effects
      if (!is.null(eqtl_effects)) {
        ct_result <- merge(ct_result, eqtl_effects, by = "variant_id", all.x = TRUE)
      } else {
        ct_result[, ":="(eqtl_genes = NA_character_, n_eqtl_genes = 0L)]
      }

      return(ct_result)
    }
    return(NULL)
  }), fill = TRUE)

  # Merge with base matrix (setkey for binary search join)
  if (nrow(qtl_effects) > 0) {
    data.table::setkeyv(base_matrix, c("variant_id", "cell_type"))
    data.table::setkeyv(qtl_effects, c("variant_id", "cell_type"))
    result <- merge(base_matrix, qtl_effects,
      by = c("variant_id", "cell_type"),
      all.x = TRUE
    )
  } else {
    result <- base_matrix
    result[, ":="(caqtl_peaks = NA_character_,
      n_caqtl_peaks = 0L,
      eqtl_genes = NA_character_,
      n_eqtl_genes = 0L)]
  }

  # Add LFSR values using pre-computed long format data
  if (!is.null(lfsr_results)) {
    # Add eQTL LFSR
    if (!is.null(lfsr_results$eqtl_long) && nrow(lfsr_results$eqtl_long) > 0) {
      eqtl_lfsr <- extract_variant_lfsr_from_long(
        result, lfsr_results$eqtl_long, "eqtl_genes"
      )
      data.table::setkeyv(result, c("variant_id", "cell_type"))
      data.table::setkeyv(eqtl_lfsr, c("variant_id", "cell_type"))
      result <- merge(result, eqtl_lfsr,
        by = c("variant_id", "cell_type"),
        all.x = TRUE
      )
    }

    # Add caQTL LFSR
    if (!is.null(lfsr_results$caqtl_long) && nrow(lfsr_results$caqtl_long) > 0) {
      caqtl_lfsr <- extract_variant_lfsr_from_long(
        result, lfsr_results$caqtl_long, "caqtl_peaks"
      )
      data.table::setkeyv(result, c("variant_id", "cell_type"))
      data.table::setkeyv(caqtl_lfsr, c("variant_id", "cell_type"))
      result <- merge(result, caqtl_lfsr,
        by = c("variant_id", "cell_type"),
        all.x = TRUE
      )
    }
  }

  # Fill missing values
  result[is.na(n_caqtl_peaks), n_caqtl_peaks := 0L]
  result[is.na(n_eqtl_genes), n_eqtl_genes := 0L]

  data.table::setorder(result, variant_id, cell_type)
  return(result)
}

#' Extract LFSR values for variants
#'
#' Helper function to extract LFSR values for variant-feature pairs
#'
#' @param variant_matrix Variant matrix with QTL effects
#' @param lfsr_data LFSR data table
#' @param cell_types Vector of cell type column names
#' @param feature_col Column containing affected features
#' @return data.table with LFSR values
#' @importFrom data.table setDT setnames
#' @keywords internal
extract_variant_lfsr <- function(variant_matrix, lfsr_data, cell_types, feature_col) {
  data.table::setDT(lfsr_data)
  data.table::setDT(variant_matrix)

  lfsr_col_name <- paste0(feature_col, "_lfsr")

  # Identify rows with non-empty features
  has_features <- variant_matrix[
    !is.na(get(feature_col)) & nchar(get(feature_col)) > 0,
    .(variant_id, cell_type, features = get(feature_col))
  ]

  if (nrow(has_features) == 0L) {
    result <- variant_matrix[, .(variant_id, cell_type)]
    result[, (lfsr_col_name) := NA_character_]
    return(result)
  }

  # Identify LFSR columns present in lfsr_data
  lfsr_cols <- intersect(cell_types, names(lfsr_data))
  if (length(lfsr_cols) == 0L) {
    result <- variant_matrix[, .(variant_id, cell_type)]
    result[, (lfsr_col_name) := NA_character_]
    return(result)
  }

  # Melt lfsr_data to long format once (replaces per-row column lookups)
  lfsr_long <- data.table::melt(lfsr_data,
    id.vars = c("feature_id", "variant_id"),
    measure.vars = lfsr_cols,
    variable.name = "cell_type",
    value.name = "lfsr_value"
  )
  lfsr_long[, cell_type := as.character(cell_type)]
  data.table::setkey(lfsr_long, variant_id, feature_id, cell_type)

  # Expand comma-separated features into individual rows
  expanded <- has_features[, .(feature_id = unlist(strsplit(features, ",", fixed = TRUE))),
    by = .(variant_id, cell_type)
  ]

  # Keyed join to look up LFSR values (replaces O(n×m) row scans)
  expanded[lfsr_long, lfsr_value := i.lfsr_value,
    on = .(variant_id, feature_id, cell_type)
  ]

  # Aggregate: paste non-NA LFSR values per (variant_id, cell_type)
  agg <- expanded[,
    {
      vals <- lfsr_value[!is.na(lfsr_value)]
      list(lfsr_str = if (length(vals) == 0L) NA_character_ else paste(round(vals, 4), collapse = ","))
    },
    by = .(variant_id, cell_type)
  ]

  # Build full result with all rows from variant_matrix
  result <- variant_matrix[, .(variant_id, cell_type)]
  result[agg, (lfsr_col_name) := i.lfsr_str, on = .(variant_id, cell_type)]

  return(result)
}

#' Extract LFSR values for variants from pre-computed long format
#'
#' @param variant_matrix Variant matrix with QTL effects
#' @param lfsr_long Pre-computed long-format LFSR data.table [feature_id, variant_id, cell_type, lfsr]
#' @param feature_col Column containing affected features (e.g., "eqtl_genes", "caqtl_peaks")
#' @return data.table with LFSR values
#' @keywords internal
extract_variant_lfsr_from_long <- function(variant_matrix, lfsr_long, feature_col) {
  data.table::setDT(variant_matrix)

  lfsr_col_name <- paste0(feature_col, "_lfsr")

  # Identify rows with non-empty features
  has_features <- variant_matrix[
    !is.na(get(feature_col)) & nchar(get(feature_col)) > 0,
    .(variant_id, cell_type, features = get(feature_col))
  ]

  if (nrow(has_features) == 0L) {
    result <- variant_matrix[, .(variant_id, cell_type)]
    result[, (lfsr_col_name) := NA_character_]
    return(result)
  }

  # Expand comma-separated features into individual rows
  expanded <- has_features[, .(feature_id = unlist(strsplit(features, ",", fixed = TRUE))),
    by = .(variant_id, cell_type)
  ]

  # Keyed join to look up LFSR values from pre-computed long format
  # lfsr_long has columns: feature_id, variant_id, cell_type, lfsr
  expanded[lfsr_long, lfsr_value := i.lfsr,
    on = .(variant_id, feature_id, cell_type)
  ]

  # Aggregate: paste non-NA LFSR values per (variant_id, cell_type)
  agg <- expanded[,
    {
      vals <- lfsr_value[!is.na(lfsr_value)]
      list(lfsr_str = if (length(vals) == 0L) NA_character_ else paste(round(vals, 4), collapse = ","))
    },
    by = .(variant_id, cell_type)
  ]

  # Build full result with all rows from variant_matrix
  result <- variant_matrix[, .(variant_id, cell_type)]
  result[agg, (lfsr_col_name) := i.lfsr_str, on = .(variant_id, cell_type)]

  return(result)
}

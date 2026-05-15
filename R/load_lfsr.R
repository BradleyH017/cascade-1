#' Load pre-computed LFSR results from external files
#'
#' @param config Configuration list containing file paths and cell types
#' @return List containing LFSR data tables for eQTL and caQTL
load_lfsr_results <- function(config) {
  cli::cli_alert_info("Loading pre-computed LFSR values")

  # Load eQTL LFSR
  lfsr_mapping <- if (!is.null(config$column_mapping$lfsr)) config$column_mapping$lfsr else NULL
  eqtl_lfsr_data <- NULL
  if (!is.null(config$file_patterns$eqtl_lfsr)) {
    validate_file_exists(config$file_patterns$eqtl_lfsr, "eQTL LFSR file")
    cli::cli_alert_info("Loading eQTL LFSR from: {.file {config$file_patterns$eqtl_lfsr}}")
    eqtl_lfsr_data <- load_and_parse_lfsr_file(config$file_patterns$eqtl_lfsr, lfsr_mapping)
    cli::cli_alert_success("Loaded {nrow(eqtl_lfsr_data)} eQTL feature-variant pairs")
  }

  # Load caQTL LFSR
  caqtl_lfsr_data <- NULL
  if (!is.null(config$file_patterns$caqtl_lfsr)) {
    validate_file_exists(config$file_patterns$caqtl_lfsr, "caQTL LFSR file")
    cli::cli_alert_info("Loading caQTL LFSR from: {.file {config$file_patterns$caqtl_lfsr}}")
    caqtl_lfsr_data <- load_and_parse_lfsr_file(config$file_patterns$caqtl_lfsr, lfsr_mapping)
    cli::cli_alert_success("Loaded {nrow(caqtl_lfsr_data)} caQTL peak-variant pairs")
  }

  # Detect cell type columns using hierarchy-based detection
  hierarchy <- if (!is.null(config$hierarchy)) config$hierarchy else DEFAULT_CELL_HIERARCHY
  eqtl_ct_cols <- if (!is.null(eqtl_lfsr_data) && nrow(eqtl_lfsr_data) > 0) {
    detect_ct_columns(eqtl_lfsr_data, hierarchy)
  } else {
    character(0)
  }
  caqtl_ct_cols <- if (!is.null(caqtl_lfsr_data) && nrow(caqtl_lfsr_data) > 0) {
    detect_ct_columns(caqtl_lfsr_data, hierarchy)
  } else {
    character(0)
  }

  # Build long-format and feature-level aggregated tables
  eqtl_long <- if (!is.null(eqtl_lfsr_data) && length(eqtl_ct_cols) > 0) {
    melt_lfsr_to_long(eqtl_lfsr_data, eqtl_ct_cols)
  } else {
    data.table::data.table(
      feature_id = character(), variant_id = character(),
      cell_type = character(), lfsr = numeric()
    )
  }

  caqtl_long <- if (!is.null(caqtl_lfsr_data) && length(caqtl_ct_cols) > 0) {
    melt_lfsr_to_long(caqtl_lfsr_data, caqtl_ct_cols)
  } else {
    data.table::data.table(
      feature_id = character(), variant_id = character(),
      cell_type = character(), lfsr = numeric()
    )
  }

  # Pass L2→L1 mapping so feature lookup keys are always L1
  hierarchy <- if (!is.null(config$hierarchy)) config$hierarchy else DEFAULT_CELL_HIERARCHY
  l2_map <- get_mapping_columns(hierarchy)
  eqtl_feature <- aggregate_lfsr_by_feature(eqtl_long, l2_to_l1 = l2_map)
  caqtl_feature <- aggregate_lfsr_by_feature(caqtl_long, l2_to_l1 = l2_map)

  return(list(
    eqtl_feature = eqtl_feature,
    caqtl_feature = caqtl_feature,
    eqtl_long = eqtl_long,
    caqtl_long = caqtl_long
  ))
}

#' Melt wide LFSR data.table to long format
#'
#' @param wide_dt Wide LFSR data.table with feature_id, variant_id, and CT columns
#' @param cell_type_cols Character vector of cell type column names
#' @return Long-format data.table [feature_id, variant_id, cell_type, lfsr]
#' @keywords internal
melt_lfsr_to_long <- function(wide_dt, cell_type_cols) {
  ct_cols_present <- intersect(cell_type_cols, names(wide_dt))
  if (length(ct_cols_present) == 0) {
    return(data.table::data.table(
      feature_id = character(), variant_id = character(),
      cell_type = character(), lfsr = numeric()
    ))
  }
  long <- data.table::melt(wide_dt,
    id.vars = c("feature_id", "variant_id"),
    measure.vars = ct_cols_present,
    variable.name = "cell_type",
    value.name = "lfsr",
    variable.factor = FALSE
  )
  # Drop NA rows to reduce memory
  long <- long[!is.na(lfsr)]
  data.table::setkey(long, feature_id, variant_id, cell_type)
  long
}

#' Aggregate long LFSR to feature-level lookup
#'
#' Maps L2 cell types to L1 before aggregation so that the resulting lookup
#' keys match the L1 names used by the C++ categorization function.
#'
#' @param long_dt Long LFSR data.table [feature_id, variant_id, cell_type, lfsr]
#' @param l2_to_l1 Named list mapping L2 → L1 cell type column names
#' @return Feature-level lookup data.table [feature_id, cell_type, min_lfsr] with L1 keys
#' @keywords internal
aggregate_lfsr_by_feature <- function(long_dt, l2_to_l1 = NULL) {
  if (nrow(long_dt) == 0) {
    return(data.table::data.table(
      feature_id = character(), cell_type = character(), min_lfsr = numeric()
    ))
  }
  # Map L2 cell types to L1 for consistent keys with C++ categorization
  dt <- data.table::copy(long_dt)
  if (!is.null(l2_to_l1) && length(l2_to_l1) > 0) {
    l2_names <- names(l2_to_l1)
    dt[cell_type %in% l2_names, cell_type := unlist(l2_to_l1[cell_type])]
  }
  # Aggregate: min LFSR per feature × L1 cell type
  feature_lookup <- dt[, .(min_lfsr = min(lfsr, na.rm = TRUE)), by = .(feature_id, cell_type)]
  feature_lookup[is.infinite(min_lfsr), min_lfsr := NA_real_]
  data.table::setkey(feature_lookup, feature_id, cell_type)
  feature_lookup
}

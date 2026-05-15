#' @title Internal Utility Functions
#' @description Shared helpers to reduce code duplication across the package.

#' Quietly run a system command and return character output
#'
#' Returns \code{character(0)} on any failure. Suppresses both R warnings
#' and the subprocess's stderr so we can probe optional system facilities
#' (e.g., \code{/proc/meminfo}) without polluting test output.
#'
#' @param cmd Shell command to run
#' @noRd
.system_quiet <- function(cmd) {
  out <- tryCatch(
    suppressWarnings(system(cmd, intern = TRUE, ignore.stderr = TRUE)),
    error = function(e) character(0)
  )
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    return(character(0))
  }
  out
}

#' Estimate the current process's resident set size in MB
#'
#' Falls back to a sensible default when the platform does not expose
#' an obvious mechanism (e.g., Windows). Uses \code{ps -o rss=} on
#' Linux and macOS.
#'
#' @param fallback_mb Numeric default to return if measurement fails
#' @return Numeric MB
#' @keywords internal
current_rss_mb <- function(fallback_mb = 2000) {
  out <- .system_quiet(paste("ps -o rss= -p", Sys.getpid()))
  if (length(out) == 0L) {
    return(fallback_mb)
  }
  val <- suppressWarnings(as.numeric(out)) / 1024
  if (length(val) == 0L || !is.finite(val)) fallback_mb else val
}

#' Estimate available system memory in MB
#'
#' Linux uses \code{/proc/meminfo}; macOS uses \code{sysctl hw.memsize}
#' (which reports total physical memory rather than currently-available,
#' which is a conservative upper bound for forking heuristics).
#'
#' @param fallback_mb Numeric default to return if measurement fails
#' @return Numeric MB
#' @keywords internal
total_mem_mb <- function(fallback_mb = 32000) {
  sysname <- Sys.info()[["sysname"]]
  if (sysname == "Linux") {
    out <- .system_quiet("grep MemAvailable /proc/meminfo")
    if (length(out) == 0L) {
      return(fallback_mb)
    }
    val <- suppressWarnings(as.numeric(gsub("[^0-9]", "", out))) / 1024
    if (length(val) == 0L || !is.finite(val)) fallback_mb else val
  } else if (sysname == "Darwin") {
    out <- .system_quiet("sysctl -n hw.memsize")
    if (length(out) == 0L) {
      return(fallback_mb)
    }
    val <- suppressWarnings(as.numeric(out)) / 1024 / 1024
    if (length(val) == 0L || !is.finite(val)) fallback_mb else val
  } else {
    fallback_mb
  }
}

#' Validate that a file exists, aborting with a descriptive message if not
#'
#' URLs are accepted as-is (data.table::fread can stream them); only on-disk
#' paths are checked with file.exists.
#'
#' @param file_path Path to check
#' @param description Human-readable description of what the file is
#' @keywords internal
validate_file_exists <- function(file_path, description) {
  if (R.utils::isUrl(file_path)) {
    return(invisible(NULL))
  }
  if (!file.exists(file_path)) {
    cli::cli_abort("{description} not found: {.file {file_path}}")
  }
}

#' Rename columns in a data.table using a mapping
#'
#' Renames input columns to internal names. NULL entries are treated as
#' optional and skipped. Fails with a clear error if a required mapped
#' column is missing.
#'
#' @param dt A data.table
#' @param mapping Named list: internal_name = "input_column_name". NULL values are optional.
#' @keywords internal
rename_columns <- function(dt, mapping) {
  for (internal_name in names(mapping)) {
    input_name <- mapping[[internal_name]]
    if (is.null(input_name)) next
    if (!input_name %in% names(dt)) {
      stop(
        "Expected column '", input_name, "' (mapped to '", internal_name,
        "') not found. Available: ",
        paste(head(names(dt), 15), collapse = ", ")
      )
    }
    if (input_name != internal_name) {
      data.table::setnames(dt, input_name, internal_name)
    }
  }
}

#' Detect cell type columns using hierarchy
#'
#' Returns column names present in the data that match cell types defined
#' in the hierarchy (L1, bulk, other, and mapped lower-level types).
#'
#' @param dt A data.table
#' @param hierarchy A CellTypeHierarchy object
#' @return Character vector of matching column names
#' @keywords internal
detect_ct_columns <- function(dt, hierarchy) {
  all_ct_cols <- unique(c(
    hierarchy$l1_columns,
    hierarchy$bulk_columns,
    hierarchy$other_columns,
    names(get_mapping_columns(hierarchy))
  ))
  intersect(all_ct_cols, names(dt))
}

#' Load and parse an LFSR file, splitting the id column into feature_id and variant_id
#'
#' @param file_path Path to the LFSR TSV file
#' @return data.table keyed by (feature_id, variant_id)
#' @importFrom data.table fread tstrsplit setkey
#' @keywords internal
load_and_parse_lfsr_file <- function(file_path, column_mapping = NULL) {
  dt <- data.table::fread(file_path, header = TRUE, sep = "\t")
  id_col <- if (!is.null(column_mapping)) column_mapping$id_column else "id"
  id_sep <- if (!is.null(column_mapping)) column_mapping$id_separator else ":"
  dt[, c("feature_id", "variant_id") := data.table::tstrsplit(get(id_col), id_sep, fixed = TRUE)]
  dt[, (id_col) := NULL]
  data.table::setkey(dt, feature_id, variant_id)
  dt
}

#' Extract QTL effects for a single cell type's data
#'
#' Aggregates variant-feature pairs, collapsing feature IDs per variant.
#' Expects standardized column names ("variant_id", "feature_id") produced
#' by \code{rename_columns()} in the loading pipeline.
#'
#' @param qtl_data data.table of QTL results for one cell type (with
#'   standardized column names: "variant_id", "feature_id")
#' @param feature_col_out Output column name for collapsed features (e.g., "caqtl_peaks")
#' @param count_col_out Output column name for feature count (e.g., "n_caqtl_peaks")
#' @return data.table with variant_id, feature list, and count; or NULL
#' @keywords internal
extract_qtl_effects <- function(qtl_data, feature_col_out, count_col_out) {
  if (is.null(qtl_data) || nrow(qtl_data) == 0L) {
    return(NULL)
  }

  # Data should already have standardized column names from rename_columns()
  if (!"variant_id" %in% names(qtl_data)) {
    stop("No 'variant_id' column found. Data should be standardized by rename_columns() before calling extract_qtl_effects().")
  }
  if (!"feature_id" %in% names(qtl_data)) {
    stop("No 'feature_id' column found. Data should be standardized by rename_columns() before calling extract_qtl_effects().")
  }

  effects <- qtl_data[, .(
    features = paste(unique(feature_id), collapse = ","),
    n_features = length(unique(feature_id))
  ), by = "variant_id"]

  data.table::setnames(
    effects, c("features", "n_features"),
    c(feature_col_out, count_col_out)
  )
  effects
}

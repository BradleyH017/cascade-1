#' @title Metadata and Peak Loading Functions
#' @description Functions for loading peak data, meta-analysis data, and CS cluster mappings

#' Load Peak Data
#'
#' Load and process peak data for overlap detection
#'
#' @param peak_file Path to ATAC peak BED file
#' @param target_chrs Target chromosomes (NULL for all chromosomes)
#' @param quiet If TRUE, suppress informational logging
#' @return List with data table and GRanges object
#' @noRd
load_peaks <- function(peak_file, target_chrs = NULL, quiet = FALSE) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
    !requireNamespace("IRanges", quietly = TRUE)) {
    cli::cli_abort("GenomicRanges and IRanges packages are required for peak overlap detection")
  }

  # Don't use progress step here as this function may be called from within loops

  # Auto-detect file format based on first few lines
  first_lines <- readLines(peak_file, n = 5)

  # Check if it's text format (chr-start-end format)
  text_format <- all(grepl("^chr[^-]+-[0-9]+-[0-9]+", first_lines))

  if (text_format) {
    if (!quiet) cli::cli_alert_info("Detected text format (chr-start-end)")

    # Use fread for loading
    peak_data <- data.table::fread(peak_file, header = FALSE, col.names = "peak_id")

    # Vectorized parsing using data.table operations
    # Split all peak IDs at once
    coords <- data.table::tstrsplit(peak_data$peak_id, "-", type.convert = TRUE)

    peaks <- data.table::data.table(
      chr = coords[[1]],
      start = as.integer(coords[[2]]) - 1, # Convert to 0-based for consistency with BED processing
      end = as.integer(coords[[3]]),
      peak_id = peak_data$peak_id,
      score = 1000, # Default score
      strand = "*" # Default strand
    )
  } else {
    if (!quiet) cli::cli_alert_info("Detected BED format")

    # Load BED format - standard BED format
    peaks <- data.table::fread(peak_file,
      header = FALSE,
      col.names = c("chr", "start", "end", "peak_id", "score", "strand")
    )
  }

  # Filter by chromosomes if specified
  if (!is.null(target_chrs)) {
    peaks <- peaks[peaks$chr %in% target_chrs, ]
  }

  # Create unique peak IDs based on coordinates
  # Convert BED 0-based start to 1-based to match caQTL format
  peaks$simple_peak_id <- paste(peaks$chr, peaks$start + 1, peaks$end, sep = "-")

  # Create GRanges object for overlap queries
  gr_peaks <- GenomicRanges::GRanges(
    seqnames = peaks$chr,
    ranges = IRanges::IRanges(start = peaks$start + 1, end = peaks$end), # BED is 0-based
    peak_id = peaks$simple_peak_id,
    original_id = peaks$peak_id
  )

  if (!quiet) cli::cli_alert_success("Loaded {.val {length(gr_peaks)}} peaks")
  return(list(dt = peaks, gr = gr_peaks))
}

#' Load Meta Data with Pre-computed Cochran's Q Values
#'
#' Load meta data containing pre-computed Cochran's Q heterogeneity p-values
#'
#' @param meta_file Path to the meta data file
#' @param cache Cache object for memoization (optional)
#' @param column_mapping Named list mapping internal names to input column names.
#'   If NULL, uses DEFAULT_COLUMN_MAPPING$meta.
#' @return Data table with variant-phenotype pairs and heterogeneity p-values
load_meta_data <- function(meta_file, cache = NULL, column_mapping = NULL) {
  validate_file_exists(meta_file, "Meta data file")

  cli::cli_alert_info("Loading meta file")

  mapping <- if (!is.null(column_mapping)) column_mapping else DEFAULT_COLUMN_MAPPING$meta

  # Load only the required columns (skips ~45% of data: meta_beta, meta_se, meta_nlog10p, meta_q, cell_types)
  select_cols <- unname(unlist(mapping))
  meta_data <- cached_fread(meta_file, cache = cache, select = select_cols)
  if (!all(select_cols %in% names(meta_data))) {
    missing <- setdiff(select_cols, names(meta_data))
    cli::cli_abort("Meta file missing required columns: {.val {missing}}. Expected columns from mapping: {.val {select_cols}}")
  }

  # Rename columns to internal names

  rename_columns(meta_data, mapping)

  # Convert nlog10p to p-value for consistency with existing code
  # cochran_q_nlog10p is -log10(p), so p = 10^(-cochran_q_nlog10p)
  meta_data$cochran_q_pval <- 10^(-meta_data$cochran_q_nlog10p)

  # Validate direction column values
  valid_directions <- c("+", "-", "+-", "-+", NA)
  invalid_dirs <- setdiff(unique(meta_data$direction), valid_directions)
  if (length(invalid_dirs) > 0) {
    cli::cli_abort("Invalid direction values found: {.val {invalid_dirs}}. Valid values are: +, -, +-, -+, NA")
  }

  # Validate max_pip column
  if (any(meta_data$max_pip > 1 | meta_data$max_pip < 0, na.rm = TRUE)) {
    cli::cli_abort("max_pip values must be between 0 and 1")
  }

  # Validate max_chisq column (chi-square values should be non-negative)
  if (any(meta_data$max_chisq < 0, na.rm = TRUE)) {
    cli::cli_abort("max_chisq values must be non-negative")
  }

  # Set multiple indices for different access patterns
  data.table::setkey(meta_data, feature_id, variant_id)

  # Create secondary indices for fast filtering
  data.table::setindex(meta_data, max_pip)
  data.table::setindex(meta_data, max_chisq)
  data.table::setindex(meta_data, feature_id)
  data.table::setindex(meta_data, variant_id)

  dir_summary <- table(meta_data$direction, useNA = "ifany")
  cli::cli_alert_info("Direction summary: {.val {paste(names(dir_summary), '=', dir_summary, collapse = ', ')}}")
  cli::cli_alert_info("Max PIP range: [{.val {round(min(meta_data$max_pip, na.rm = TRUE), 3)}}, {.val {round(max(meta_data$max_pip, na.rm = TRUE), 3)}}]")
  cli::cli_alert_info("Max Chi-square range: [{.val {round(min(meta_data$max_chisq, na.rm = TRUE), 2)}}, {.val {round(max(meta_data$max_chisq, na.rm = TRUE), 2)}}]")

  cli::cli_alert_success("Loaded {.val {nrow(meta_data)}} variant-phenotype pairs")

  return(meta_data)
}

#' Load CS cluster mapping data
#'
#' Loads the CS cluster file that maps credible sets to clusters across cell types
#'
#' @param cs_cluster_file Path to the CS cluster file
#' @param cache Cache object for memoization (optional)
#' @param column_mapping Named list mapping internal names to input column names.
#'   If NULL, uses DEFAULT_COLUMN_MAPPING$cs_clusters.
#' @return Data table with CS to cluster mappings indexed by feature
load_cs_clusters <- function(cs_cluster_file, cache = NULL, column_mapping = NULL) {
  if (!file.exists(cs_cluster_file)) {
    cli::cli_alert_warning("CS cluster file not found: {.file {cs_cluster_file}}")
    return(NULL)
  }

  cli::cli_alert_info("Loading CS cluster mappings")

  mapping <- if (!is.null(column_mapping)) column_mapping else DEFAULT_COLUMN_MAPPING$cs_clusters

  # Load the CS cluster data
  cs_clusters <- cached_fread(cs_cluster_file, cache = cache)

  # Check for required columns (input names from mapping)
  required_cols <- unlist(mapping)
  if (!all(required_cols %in% names(cs_clusters))) {
    missing <- setdiff(required_cols, names(cs_clusters))
    cli::cli_abort("CS cluster file missing required columns: {.val {missing}}")
  }

  # Rename columns to internal names
  rename_columns(cs_clusters, mapping)

  # Set multiple keys for efficient access
  data.table::setkey(cs_clusters, feature_id, cluster_id)
  data.table::setindex(cs_clusters, cell_type)
  data.table::setindex(cs_clusters, qtl_type)

  # Print summary
  n_features <- length(unique(cs_clusters$feature_id))
  n_clusters <- length(unique(cs_clusters$cluster_id))
  n_cell_types <- length(unique(cs_clusters$cell_type))

  cli::cli_alert_success("Loaded CS clusters: {.val {n_features}} features, {.val {n_clusters}} clusters, {.val {n_cell_types}} cell types")

  return(cs_clusters)
}

#' Load CS cluster variant data
#'
#' Loads the file mapping clusters to their constituent variants
#'
#' @param cs_cluster_variant_file Path to the CS cluster variant file
#' @param cache Cache object for memoization (optional)
#' @param column_mapping Named list mapping internal names to input column names.
#'   If NULL, uses DEFAULT_COLUMN_MAPPING$cs_cluster_variants.
#' @return Data table with cluster to variant mappings
load_cs_cluster_variants <- function(cs_cluster_variant_file, cache = NULL, column_mapping = NULL) {
  if (!file.exists(cs_cluster_variant_file)) {
    cli::cli_alert_warning("CS cluster variant file not found: {.file {cs_cluster_variant_file}}")
    return(NULL)
  }

  cli::cli_alert_info("Loading CS cluster variant mappings")

  mapping <- if (!is.null(column_mapping)) column_mapping else DEFAULT_COLUMN_MAPPING$cs_cluster_variants

  # Load the CS cluster variant data
  cs_cluster_variants <- cached_fread(cs_cluster_variant_file, cache = cache)

  # Check for required columns (input names from mapping)
  required_cols <- unlist(mapping)
  if (!all(required_cols %in% names(cs_cluster_variants))) {
    missing <- setdiff(required_cols, names(cs_cluster_variants))
    cli::cli_abort("CS cluster variant file missing required columns: {.val {missing}}")
  }

  # Rename columns to internal names
  rename_columns(cs_cluster_variants, mapping)

  # Ensure data.table format and set indices
  if (!data.table::is.data.table(cs_cluster_variants)) {
    data.table::setDT(cs_cluster_variants)
  }

  # Set keys for efficient access
  data.table::setkey(cs_cluster_variants, cluster_id, variant_id)
  data.table::setindex(cs_cluster_variants, feature_ids)

  # Print summary
  n_clusters <- length(unique(cs_cluster_variants$cluster_id))
  n_variants <- length(unique(cs_cluster_variants$variant_id))
  n_features <- length(unique(cs_cluster_variants$feature_ids))

  cli::cli_alert_success("Loaded cluster variants: {.val {n_clusters}} clusters, {.val {n_variants}} variants, {.val {n_features}} features")

  return(cs_cluster_variants)
}

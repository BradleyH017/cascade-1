#' @title Feature Data Loading Functions
#' @description Functions for loading various data types required for categorization

#' Construct File Paths
#'
#' Helper function to build cell type and chromosome-specific file paths
#' by replacing placeholders in pattern strings
#'
#' @param pattern File pattern with placeholders {CELL_TYPE} and/or {CHR}
#' @param cell_type Cell type to substitute for {CELL_TYPE}
#' @param chromosome Chromosome to substitute for {CHR} (optional)
#' @return Modified file path with placeholders replaced
#' @keywords internal
construct_file_paths <- function(pattern, cell_type = NULL, chromosome = NULL) {
  path <- pattern

  # Replace cell type placeholder if provided
  if (!is.null(cell_type)) {
    path <- gsub("{CELL_TYPE}", cell_type, path, fixed = TRUE)
  }

  # Replace chromosome placeholder if provided
  if (!is.null(chromosome)) {
    path <- gsub("{CHR}", chromosome, path, fixed = TRUE)
  }

  return(path)
}

#' Load Feature Data
#'
#' Generic function to load feature data (genes or peaks) with ACAT results
#'
#' @param config Configuration list with file patterns
#' @param feature_type Type of feature ("gene" or "peak")
#' @param chromosomes Chromosomes to analyze
#' @param num_cores Number of cores for parallel processing (NULL to use config$parameters$n_cores)
#' @return List with feature data
load_feature_data <- function(config, feature_type, chromosomes = NULL, num_cores = NULL) {
  cli::cli_alert_info("Loading {feature_type} data")

  # Check config type
  if (!is.list(config)) {
    stop("Config must be a list object, not ", class(config))
  }

  # Initialize cache
  cache <- init_cache(config)

  # Extract patterns and cell types based on feature type
  if (feature_type == "gene") {
    acat_pattern <- config$file_patterns$eqtl_acat
    susie_pattern <- config$file_patterns$eqtl_susie
    susie_col_mapping <- config$column_mapping$eqtl_susie
  } else if (feature_type == "peak") {
    acat_pattern <- config$file_patterns$caqtl_acat
    susie_pattern <- config$file_patterns$caqtl_susie
    susie_col_mapping <- config$column_mapping$caqtl_susie
  } else {
    cli::cli_abort("Unknown feature type: {feature_type}. Must be 'gene' or 'peak'.")
  }

  if (is.null(acat_pattern)) {
    cli::cli_abort("{feature_type}_acat file pattern not found in config")
  }

  if (is.null(susie_pattern)) {
    cli::cli_abort("{feature_type}_susie file pattern not found in config. SuSiE results are required for variant analysis.")
  }

  if (is.null(chromosomes)) {
    cli::cli_abort("Chromosomes must be specified in config. Cannot proceed without chromosome information.")
  }

  cell_types <- as.character(unlist(config$cell_types))

  # Determine column mapping once before the loop
  if (feature_type == "gene") {
    col_mapping_key <- "eqtl_acat"
    default_col_mapping <- list(feature_id = "phenotype_id", q_value = "ACAT_q")
  } else {
    col_mapping_key <- "caqtl_acat"
    default_col_mapping <- list(feature_id = "phenotype_id", q_value = "qval")
  }
  col_mapping <- config$column_mapping[[col_mapping_key]]
  if (is.null(col_mapping)) col_mapping <- default_col_mapping

  # Load ACAT results (only the 2 needed columns via select)
  acat_results <- list()
  for (ct in cell_types) {
    acat_file <- construct_file_paths(acat_pattern, cell_type = ct)
    if (file.exists(acat_file)) {
      cli::cli_alert_info("Loading ACAT results for {ct}")
      acat_data <- cached_fread(acat_file,
        cache = cache,
        select = c(col_mapping$feature_id, col_mapping$q_value)
      )

      # Note: Chromosome filtering is handled via SuSiE results, not at ACAT level

      # Extract data using column mapping
      acat_results[[ct]] <- data.frame(
        feature_id = acat_data[[col_mapping$feature_id]],
        q_value = acat_data[[col_mapping$q_value]],
        stringsAsFactors = FALSE
      )
    } else {
      cli::cli_abort("ACAT file not found: {.file {acat_file}}")
    }
  }

  if (length(acat_results) == 0) {
    cli::cli_abort("No ACAT files could be loaded")
  }

  # Organize by feature - get the list of features before loading SuSiE data
  all_features <- unique(unlist(lapply(acat_results, function(x) x$feature_id)))

  cli::cli_alert_info("Found {.val {length(all_features)}} {feature_type}s across all cell types")

  # Apply debug limit if enabled
  if (!is.null(config$debug) && config$debug$enabled) {
    debug_limit_name <- paste0("n_", feature_type, "s")
    n_features_limit <- config$debug[[debug_limit_name]]
    if (!is.null(n_features_limit) && length(all_features) > n_features_limit) {
      cli::cli_alert_warning("Limiting to {.val {n_features_limit}} {feature_type}s out of {.val {length(all_features)}}")
      all_features <- all_features[1:n_features_limit]
    }
  }

  # Load SuSiE results with lazy loading for identified features
  cli::cli_alert_info("Loading SuSiE results")
  # Use num_cores from parameter or config
  cores_to_use <- if (!is.null(num_cores)) num_cores else config$parameters$n_cores
  susie_results <- load_susie_data(
    susie_pattern = susie_pattern,
    cell_types = cell_types,
    chromosomes = chromosomes,
    feature_ids_needed = all_features, # Required for lazy loading
    parallel = TRUE,
    num_cores = cores_to_use,
    column_mapping = susie_col_mapping
  )

  if (is.null(susie_results) || length(susie_results) == 0) {
    cli::cli_abort("Failed to load SuSiE results. Please check that SuSiE files exist for the specified pattern and chromosomes.")
  }

  # Create a matrix of q-values for lookup
  feature_qval_matrix <- matrix(NA,
    nrow = length(all_features), ncol = length(acat_results),
    dimnames = list(all_features, names(acat_results))
  )

  # Fill the matrix using vectorized operations
  for (ct in names(acat_results)) {
    ct_data <- acat_results[[ct]]
    # Match all features at once
    match_idx <- match(all_features, ct_data$feature_id)
    valid_idx <- !is.na(match_idx)
    feature_qval_matrix[valid_idx, ct] <- ct_data$q_value[match_idx[valid_idx]]
  }

  # Convert to data frame format expected by downstream functions
  acat_df <- as.data.frame(feature_qval_matrix)
  acat_df$feature_id <- all_features

  # Ensure we only select columns that exist
  available_cols <- intersect(cell_types, colnames(acat_df))
  acat_df <- acat_df[, c("feature_id", available_cols)]

  # Create result list
  result <- list(
    acat_matrix = acat_df,
    feature_ids = all_features,
    susie_results = susie_results,
    cell_types = cell_types,
    n_features = length(all_features),
    cache = cache # Include cache object for downstream use
  )

  return(result)
}

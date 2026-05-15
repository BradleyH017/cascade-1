# Cache management for CASCADE data loading.
# Internal functions; uses memoise for filesystem-backed memoization.

#' @import memoise
#' @import digest
NULL

# Null-coalescing operator (base R does not provide %||% before R 4.4).
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Initialize Cache System
#'
#' Creates a filesystem cache using memoise for the CASCADE pipeline
#'
#' @param config Configuration list with cache settings
#' @return memoise cache object or NULL if disabled
#' @keywords internal
init_cache <- function(config) {
  if (!is.list(config)) {
    stop("Config must be a list object in init_cache, not ", class(config))
  }

  # Try to access cache safely
  cache_config <- tryCatch(
    {
      config$cache
    },
    error = function(e) {
      cli::cli_alert_danger("Error accessing config$cache: {e$message}")
      stop("$ operator is invalid for atomic vectors")
    }
  )

  # Check if cache is enabled
  if (is.null(cache_config) ||
    is.null(cache_config$enabled) ||
    !cache_config$enabled) {
    return(NULL)
  }

  # Ensure output directory exists
  if (is.null(config$output_dir)) {
    cli::cli_warn("Cache enabled but output_dir not specified in config")
    return(NULL)
  }

  output_dir <- config$output_dir

  # Create cache directory
  cache_dir <- file.path(output_dir, ".cascade_cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Create filesystem cache object with directory path attached
  cache <- memoise::cache_filesystem(cache_dir)
  attr(cache, "cache_dir") <- cache_dir

  # Check if configuration has changed
  config_hash <- generate_config_hash(config)
  hash_file <- file.path(cache_dir, ".config_hash")

  if (file.exists(hash_file)) {
    old_hash <- readLines(hash_file, n = 1, warn = FALSE)
    if (old_hash != config_hash) {
      cli::cli_alert_info("Configuration changed, clearing cache")
      clear_cache_internal(cache_dir)
      writeLines(config_hash, hash_file)
    }
  } else {
    writeLines(config_hash, hash_file)
  }

  return(cache)
}

#' Generate Configuration Hash
#'
#' Create hash of configuration elements that affect data loading
#'
#' @param config Configuration list
#' @return Character hash string
#' @keywords internal
generate_config_hash <- function(config) {
  if (!is.list(config)) {
    stop("Config must be a list object in generate_config_hash, not ", class(config))
  }

  # Extract relevant config elements for cache key
  # Include both parameters and thresholds to ensure cache invalidation
  # when any analysis parameter changes
  cache_key_elements <- list(
    file_patterns = config$file_patterns,
    cell_types = config$cell_types,
    parameters = config$parameters,
    debug = config$debug,
    cascade_version = as.character(packageVersion("cascade"))
  )

  # Generate hash
  digest::digest(cache_key_elements, algo = "sha256")
}

#' Clear Cache Internal
#'
#' Remove all cache files in directory
#'
#' @param cache_dir Cache directory
#' @keywords internal
clear_cache_internal <- function(cache_dir) {
  if (dir.exists(cache_dir)) {
    # Remove all files except config hash
    cache_files <- list.files(cache_dir, full.names = TRUE)
    cache_files <- cache_files[!grepl("\\.config_hash$", cache_files)]
    if (length(cache_files) > 0) {
      unlink(cache_files, recursive = TRUE)
    }
  }
}


#' Cached fread
#'
#' Read file with qs2-based binary caching. On first read, parses the file
#' with fread and saves the result as a .qs binary file in the cache directory.
#' On subsequent reads, loads from the binary cache (5-10x faster).
#'
#' Cache invalidation is based on the source file's modification time.
#'
#' @param file_path File to read
#' @param cache Cache object from init_cache (filesystem cache with path)
#' @param ... Additional arguments to fread
#' @return Data table
#' @keywords internal
cached_fread <- function(file_path, cache = NULL, ...) {
  if (is.null(cache)) {
    return(data.table::fread(file_path, ...))
  }

  validate_file_exists(file_path, "File")

  # Build a cache key from file path, modification time, and extra args (e.g., select)
  # Including ... args ensures different column selections produce different cache keys
  mtime <- file.info(file_path)$mtime
  extra_args <- list(...)
  cache_key <- digest::digest(list(file_path, as.numeric(mtime), extra_args), algo = "xxhash64")

  # Get cache directory from attribute set by init_cache
  cache_dir <- attr(cache, "cache_dir")
  if (is.null(cache_dir) || !dir.exists(cache_dir)) {
    return(data.table::fread(file_path, ...))
  }

  # qs2 binary cache path
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs_cache_file <- file.path(cache_dir, paste0(cache_key, ".qs"))
    if (file.exists(qs_cache_file)) {
      return(qs2::qs_read(qs_cache_file))
    }
    # Read, cache, and return
    result <- data.table::fread(file_path, ...)
    tryCatch(
      qs2::qs_save(result, qs_cache_file),
      error = function(e) NULL # Silently skip cache write on error
    )
    return(result)
  }

  # Fallback: RDS cache if qs2 not available
  rds_cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
  if (file.exists(rds_cache_file)) {
    return(readRDS(rds_cache_file))
  }
  result <- data.table::fread(file_path, ...)
  tryCatch(
    saveRDS(result, rds_cache_file),
    error = function(e) NULL
  )
  return(result)
}

#' Cached readRDS
#'
#' Read RDS file with qs2-based binary caching. On first read, loads the RDS
#' and saves as .qs for faster subsequent loads.
#'
#' @param file_path RDS file to read
#' @param cache Cache object from init_cache
#' @return The object stored in the RDS file
#' @keywords internal
cached_readRDS <- function(file_path, cache = NULL) {
  if (is.null(cache)) {
    return(readRDS(file_path))
  }

  validate_file_exists(file_path, "File")

  mtime <- file.info(file_path)$mtime
  cache_key <- digest::digest(list(file_path, as.numeric(mtime)), algo = "xxhash64")

  cache_dir <- attr(cache, "cache_dir")
  if (is.null(cache_dir) || !dir.exists(cache_dir)) {
    return(readRDS(file_path))
  }

  if (requireNamespace("qs2", quietly = TRUE)) {
    qs_cache_file <- file.path(cache_dir, paste0(cache_key, ".qs"))
    if (file.exists(qs_cache_file)) {
      return(qs2::qs_read(qs_cache_file))
    }
    result <- readRDS(file_path)
    tryCatch(
      qs2::qs_save(result, qs_cache_file),
      error = function(e) NULL
    )
    return(result)
  }

  # No qs2: just read directly (RDS -> RDS cache would be pointless)
  return(readRDS(file_path))
}

#' Clear CASCADE Cache
#'
#' Remove all cached data for a configuration
#'
#' @param config Configuration list
cascade_cache_clear <- function(config) {
  cache_dir <- file.path(config$output_dir, ".cascade_cache")

  if (dir.exists(cache_dir)) {
    clear_cache_internal(cache_dir)
    cli::cli_alert_success("Cache cleared")
  } else {
    cli::cli_alert_info("No cache found")
  }
}

#' Get CASCADE Cache Information
#'
#' Display cache statistics and information
#'
#' @param config Configuration list
cascade_cache_info <- function(config) {
  cache_dir <- file.path(config$output_dir, ".cascade_cache")

  if (!dir.exists(cache_dir)) {
    cli::cli_alert_info("No cache found")
    return(invisible(NULL))
  }

  cli::cli_h2("CASCADE Cache Information")
  cli::cli_alert_info("Location: {.path {cache_dir}}")

  # Count cache files
  cache_files <- list.files(cache_dir, recursive = TRUE)
  cache_files <- cache_files[!grepl("\\.config_hash$", cache_files)]

  if (length(cache_files) > 0) {
    # Calculate total size
    sizes <- file.info(file.path(cache_dir, cache_files))$size
    total_size <- sum(sizes, na.rm = TRUE)

    cli::cli_alert_info("Entries: {.val {length(cache_files)}}")
    cli::cli_alert_info("Total size: {.val {format(total_size, units = 'auto', standard = 'SI')}}")
  } else {
    cli::cli_alert_info("Cache is empty")
  }

  # Check config hash
  hash_file <- file.path(cache_dir, ".config_hash")
  if (file.exists(hash_file)) {
    config_hash <- readLines(hash_file, n = 1, warn = FALSE)
    cli::cli_alert_info("Config hash: {substr(config_hash, 1, 8)}...")
  }

  return(invisible(list(
    location = cache_dir,
    entries = length(cache_files),
    size = if (length(cache_files) > 0) total_size else 0
  )))
}

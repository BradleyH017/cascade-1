#' @title SuSiE Data Loading Functions
#' @description Functions for loading and processing SuSiE fine-mapping results

#' SuSiE Loading Function with Lazy Loading
#'
#' Efficiently loads SuSiE data for specific features using lazy loading
#'
#' @param susie_pattern File pattern for SuSiE results
#' @param cell_types Vector of cell type names
#' @param chromosomes Vector of chromosome names
#' @param feature_ids_needed Vector of feature IDs to load
#' @param parallel Whether to use parallel processing (default: TRUE)
#' @param num_cores Number of cores for parallel processing (default: auto-detect)
#' @param column_mapping Named list mapping internal names to input column names
#'   (e.g., \code{list(feature_id = "molecular_trait_id", pip = "prob")}).
#'   If NULL, uses \code{DEFAULT_COLUMN_MAPPING$eqtl_susie}.
#'
#' @return List of SuSiE results organized by cell type and gene/peak
#' @keywords internal
load_susie_data <- function(susie_pattern, cell_types, chromosomes,
                            feature_ids_needed, parallel = TRUE, num_cores = NULL,
                            column_mapping = NULL) {
  # Validate inputs
  if (is.null(feature_ids_needed) || length(feature_ids_needed) == 0) {
    cli::cli_abort("feature_ids_needed must be provided for SuSiE data loading")
  }

  cli::cli_alert_info("Loading SuSiE data for {.val {length(feature_ids_needed)}} features")

  # Memory-aware core limit for mclapply (forking inherits full RSS).
  if (is.null(num_cores)) {
    per_worker_mb <- current_rss_mb() + 1000 # RSS + per-chromosome overhead
    safe_cores <- max(1L, floor(total_mem_mb() * 0.8 / per_worker_mb))
    num_cores <- min(safe_cores, length(chromosomes), 4L)
  }

  # Decide whether to use parallel processing
  use_parallel <- parallel && requireNamespace("parallel", quietly = TRUE) &&
    num_cores > 1 && (length(chromosomes) > 1 || length(cell_types) > 1)

  # Process chromosomes either in parallel or sequential
  if (use_parallel) {
    cli::cli_alert_info("Using parallel processing with {.val {num_cores}} cores")

    all_results <- parallel::mclapply(chromosomes, function(chr) {
      load_susie_chromosome(susie_pattern, cell_types, chr, feature_ids_needed,
        column_mapping = column_mapping
      )
    }, mc.cores = num_cores, mc.preschedule = FALSE)
  } else {
    cli::cli_alert_info("Using sequential processing")

    all_results <- lapply(chromosomes, function(chr) {
      load_susie_chromosome(susie_pattern, cell_types, chr, feature_ids_needed,
        column_mapping = column_mapping
      )
    })
  }

  # Combine results across chromosomes (same for both parallel and sequential)
  susie_results <- combine_chromosome_results(all_results, cell_types)

  return(susie_results)
}

#' Process SuSiE Data with Lazy Loading
#'
#' Load and filter SuSiE data for specific features only (lazy loading)
#'
#' @param susie_file Path to SuSiE file
#' @param feature_ids_needed Vector of feature IDs to filter
#' @param column_mapping Named list mapping internal names to input column names.
#'   If NULL, uses \code{DEFAULT_COLUMN_MAPPING$eqtl_susie}.
#'
#' @return Filtered data.table with columns renamed to internal names
#' @keywords internal
process_susie_data_lazy <- function(susie_file, feature_ids_needed,
                                    column_mapping = NULL) {
  # Determine column mapping
  if (is.null(column_mapping)) {
    column_mapping <- DEFAULT_COLUMN_MAPPING$eqtl_susie
  }

  # Build reverse lookup: internal_name -> input_name
  # to know which input columns to select from file
  input_names <- unlist(column_mapping)

  # Read just the header to check column names
  # Handle .bgz files properly
  if (grepl("\\.bgz$", susie_file)) {
    header_line <- system(paste("zcat", shQuote(susie_file), "| head -1"), intern = TRUE)
    if (length(header_line) == 0 || nchar(header_line) == 0) {
      return(data.table::data.table())
    }
    header <- strsplit(header_line, "\t")[[1]]
  } else {
    header <- names(data.table::fread(susie_file, nrows = 0))
  }

  # Select only columns that exist in the file from the mapping
  essential_cols <- intersect(input_names, header)

  if (length(essential_cols) == 0) {
    return(data.table::data.table())
  }

  # Identify the feature_id input column from the mapping
  feature_id_input <- column_mapping$feature_id
  if (is.null(feature_id_input) || !feature_id_input %in% header) {
    cli::cli_abort("Feature ID column '{feature_id_input}' not found in {.file {susie_file}}")
  }

  # Read data, selecting only mapped columns that exist
  if (grepl("\\.bgz$", susie_file)) {
    line_count <- as.numeric(system(paste("zcat", shQuote(susie_file), "| wc -l"), intern = TRUE))
    if (line_count <= 1) {
      return(data.table::data.table())
    }

    susie_data <- data.table::fread(
      cmd = paste("zcat", shQuote(susie_file)),
      select = essential_cols,
      key = feature_id_input
    )
  } else {
    susie_data <- data.table::fread(
      susie_file,
      select = essential_cols,
      key = feature_id_input
    )
  }

  # Filter to only the requested features (using input column name, before rename)
  susie_data <- susie_data[get(feature_id_input) %in% feature_ids_needed]

  # Rename columns to internal names (modifies in place)
  # Only rename columns that were actually loaded
  loaded_mapping <- column_mapping[vapply(column_mapping, function(x) {
    !is.null(x) && x %in% names(susie_data)
  }, logical(1))]
  rename_columns(susie_data, loaded_mapping)

  return(susie_data)
}

#' Load SuSiE Data for Single Chromosome with Lazy Loading
#'
#' @param susie_pattern File pattern
#' @param cell_types Cell type vector
#' @param chr Chromosome name
#' @param feature_ids_needed Vector of feature IDs to filter
#' @param column_mapping Named list mapping internal names to input column names.
#'   If NULL, uses \code{DEFAULT_COLUMN_MAPPING$eqtl_susie}.
#'
#' @return List of SuSiE results for this chromosome
#' @keywords internal
load_susie_chromosome <- function(susie_pattern, cell_types, chr, feature_ids_needed,
                                  column_mapping = NULL) {
  chr_results <- list()

  for (ct in cell_types) {
    file_path <- construct_file_paths(susie_pattern, cell_type = ct, chromosome = chr)

    if (file.exists(file_path)) {
      # Always use lazy loading for efficiency
      # process_susie_data_lazy already renames columns to internal names
      susie_data <- process_susie_data_lazy(file_path, feature_ids_needed,
        column_mapping = column_mapping
      )

      if (nrow(susie_data) > 0) {
        # Data already has internal column names from lazy loader
        processed_data <- process_susie_data_vectorized(susie_data,
          column_mapping = column_mapping
        )
        chr_results[[ct]] <- processed_data
      }
    }
  }

  return(chr_results)
}

#' Vectorized SuSiE Data Processing
#'
#' Process SuSiE data using vectorized operations
#'
#' @param susie_data Raw SuSiE data.table (columns may already be renamed to
#'   internal names if called via \code{process_susie_data_lazy}).
#' @param column_mapping Named list mapping internal names to input column names.
#'   If NULL and columns have not been renamed, uses default internal names.
#'
#' @return List of processed SuSiE results by gene/peak
#' @keywords internal
process_susie_data_vectorized <- function(susie_data, column_mapping = NULL) {
  # If column_mapping is provided but data hasn't been renamed yet (direct call),
  # apply the rename. If data was already renamed by process_susie_data_lazy,
  # this is a no-op since columns already have internal names.
  if (!is.null(column_mapping)) {
    loaded_mapping <- column_mapping[vapply(column_mapping, function(x) {
      !is.null(x) && x %in% names(susie_data)
    }, logical(1))]
    if (length(loaded_mapping) > 0) {
      rename_columns(susie_data, loaded_mapping)
    }
  }

  # Use internal column names (set after rename_columns)
  gene_col <- "feature_id"
  cs_col <- "cs_id"
  variant_col <- "variant_id"
  pip_col <- "pip"
  beta_col <- if ("beta" %in% names(susie_data)) "beta" else NULL
  se_col <- if ("se" %in% names(susie_data)) "se" else NULL

  # Truly vectorized processing using data.table for speed
  if (!data.table::is.data.table(susie_data)) {
    data.table::setDT(susie_data)
  }

  # Filter to only rows with valid credible sets to reduce processing
  valid_rows <- !is.na(susie_data[[cs_col]]) & susie_data[[cs_col]] > 0
  if (sum(valid_rows) == 0) {
    return(list()) # No valid credible sets
  }

  # Split by gene using data.table
  gene_list <- split(susie_data, by = gene_col, keep.by = FALSE)

  # Process all genes using lapply (vectorized)
  gene_results <- lapply(names(gene_list), function(gene_id) {
    gene_dt <- gene_list[[gene_id]]

    if (nrow(gene_dt) == 0) {
      return(NULL)
    }

    # Extract basic data
    pips <- gene_dt[[pip_col]]
    variants <- gene_dt[[variant_col]]
    cs_values <- gene_dt[[cs_col]]
    betas <- if (!is.null(beta_col)) gene_dt[[beta_col]] else NULL
    ses <- if (!is.null(se_col)) gene_dt[[se_col]] else NULL

    # Process credible sets vectorized
    cs_list <- NULL
    valid_cs_mask <- !is.na(cs_values) & cs_values > 0

    if (any(valid_cs_mask)) {
      valid_cs <- cs_values[valid_cs_mask]
      cs_numbers <- unique(valid_cs)

      # Create credible set indices vectorized
      cs_list <- lapply(cs_numbers, function(cs_num) {
        which(cs_values == cs_num)
      })
      names(cs_list) <- paste0("L", cs_numbers)
    }

    result <- list(
      pip = pips,
      sets = list(cs = cs_list),
      variant_names = variants
    )

    if (!is.null(betas)) {
      result$beta <- betas
    }

    if (!is.null(ses)) {
      result$se <- ses
    }

    return(result)
  })

  # Set names
  names(gene_results) <- names(gene_list)

  # Remove NULL entries
  gene_results <- gene_results[!sapply(gene_results, is.null)]

  return(gene_results)
}

#' Combine SuSiE Results from Multiple Chromosomes
#'
#' Combines chromosome results from either parallel or sequential processing
#'
#' @param chr_results_list List of chromosome results
#' @param cell_types Cell type names
#'
#' @return Combined results structure with all chromosomes merged
#' @keywords internal
combine_chromosome_results <- function(chr_results_list, cell_types) {
  # Initialize result structure for each cell type
  susie_results <- lapply(cell_types, function(ct) list())
  names(susie_results) <- cell_types

  # Combine all chromosome results
  for (chr_result in chr_results_list) {
    if (!is.null(chr_result)) {
      for (ct in names(chr_result)) {
        if (ct %in% cell_types) {
          # Append chromosome's results for this cell type
          susie_results[[ct]] <- c(susie_results[[ct]], chr_result[[ct]])
        }
      }
    }
  }

  return(susie_results)
}

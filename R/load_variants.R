#' @title Variant Data Loading Functions
#' @description Functions for loading variant data with QTL information

#' Load Variant Data for All Chromosomes with Parallelization
#'
#' Load variant data across all specified chromosomes using parallel processing
#'
#' @param config Configuration list
#' @param chromosomes Chromosomes to analyze
#' @param pip_threshold Maximum PIP threshold across cell types for filtering (default: 0.5)
#' @param min_pip_threshold Minimum PIP threshold per cell type (default: 0.1)
#' @param acat_fdr_threshold FDR threshold for ACAT filtering
#' @param peak_bed_file Path to peak BED file
#' @param column_mapping Column name mapping
#' @param num_cores Number of cores to use (NULL for auto-detect)
#' @return Combined variant data from all chromosomes
load_variant_data <- function(config, chromosomes = NULL,
                              pip_threshold = 0.5,
                              min_pip_threshold = 0.1,
                              acat_fdr_threshold = 0.05,
                              peak_bed_file = NULL,
                              column_mapping = NULL,
                              num_cores = NULL) {
  # Initialize cache
  cache <- init_cache(config)

  if (is.null(chromosomes) || length(chromosomes) == 0) {
    cli::cli_abort("No chromosomes specified for variant loading")
  }

  # Determine number of cores with memory-aware limits.
  # mclapply forks the current process, so each worker inherits the full RSS.
  # We limit cores based on current RSS and available system memory.
  rss_mb <- current_rss_mb()
  mem_mb <- total_mem_mb()

  # Each forked worker inherits rss_mb, plus needs ~2GB for chromosome loading.
  per_worker_mb <- rss_mb + 2000
  safe_cores <- max(1L, floor(mem_mb * 0.8 / per_worker_mb))
  safe_cores <- min(safe_cores, length(chromosomes), 4L) # cap at 4 to be safe

  if (!is.null(num_cores)) {
    num_cores <- min(as.integer(num_cores), safe_cores)
  } else if (!is.null(config$parameters$n_cores)) {
    num_cores <- min(as.integer(config$parameters$n_cores), safe_cores)
  } else {
    num_cores <- safe_cores
  }

  cli::cli_alert_info("Loading variant data for {.val {length(chromosomes)}} chromosomes using {.val {num_cores}} cores")
  cli::cli_alert_info("Memory: current RSS {.val {round(rss_mb)}} MB, available {.val {round(mem_mb)}} MB, safe cores: {.val {safe_cores}}")

  # Pre-load chromosome-independent data ONCE (peak-gene links).
  # These files don't have {CHR} placeholders -- loading them per-chromosome
  # would re-read the same files for every chromosome.
  cli::cli_alert_info("Pre-loading peak-gene links (once for all chromosomes)")
  preloaded_peak_gene_links <- list()
  peak_gene_files <- list()
  for (ct in config$cell_types) {
    pg_file <- construct_file_paths(config$file_patterns$peak_gene_links, cell_type = ct)
    peak_gene_files[[ct]] <- pg_file
    if (file.exists(pg_file)) {
      preloaded_peak_gene_links[[ct]] <- data.table::fread(pg_file)
    }
  }

  # Helper to set up per-chromosome SuSiE file paths
  setup_chr_files <- function(chr) {
    caqtl_files <- list()
    eqtl_files <- list()
    for (ct in config$cell_types) {
      caqtl_files[[ct]] <- construct_file_paths(config$file_patterns$caqtl_susie,
        cell_type = ct, chromosome = chr
      )
      eqtl_files[[ct]] <- construct_file_paths(config$file_patterns$eqtl_susie,
        cell_type = ct, chromosome = chr
      )
    }
    list(caqtl = caqtl_files, eqtl = eqtl_files)
  }

  # Use parallel processing to load data for each chromosome
  if (num_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    cli::cli_alert_info("Using parallel processing")
    options(mc.cores = num_cores)

    chr_results <- parallel::mclapply(chromosomes, function(chr) {
      files <- setup_chr_files(chr)
      load_variant_data_by_qtl_type(
        caqtl_files = files$caqtl,
        eqtl_files = files$eqtl,
        peak_gene_files = peak_gene_files,
        chromosomes = chr,
        pip_threshold = pip_threshold,
        min_pip_threshold = min_pip_threshold,
        acat_fdr_threshold = acat_fdr_threshold,
        peak_bed_file = config$file_patterns$peak_bed,
        column_mapping = column_mapping,
        config = config,
        quiet = TRUE,
        preloaded_peak_gene_data = preloaded_peak_gene_links
      )
    }, mc.cores = num_cores, mc.preschedule = FALSE)
  } else {
    cli::cli_alert_info("Using sequential processing")
    chr_results <- list()
    for (chr in chromosomes) {
      cli::cli_alert_info("Processing chromosome {chr}")
      files <- setup_chr_files(chr)
      chr_results[[chr]] <- load_variant_data_by_qtl_type(
        caqtl_files = files$caqtl,
        eqtl_files = files$eqtl,
        peak_gene_files = peak_gene_files,
        chromosomes = chr,
        pip_threshold = pip_threshold,
        min_pip_threshold = min_pip_threshold,
        acat_fdr_threshold = acat_fdr_threshold,
        peak_bed_file = config$file_patterns$peak_bed,
        column_mapping = column_mapping,
        config = config,
        quiet = TRUE,
        preloaded_peak_gene_data = preloaded_peak_gene_links
      )
    }
  }

  # Combine results from all chromosomes
  cli::cli_alert_info("Combining results from all chromosomes")
  combined_data <- combine_chromosome_variant_data(chr_results)

  cli::cli_alert_success("Loaded {.val {combined_data$n_variants}} variants across {.val {length(chromosomes)}} chromosomes")

  return(combined_data)
}

#' Combine Variant Data Across Chromosomes
#'
#' @param chr_results List of variant data by chromosome
#' @return Combined variant data
#' @keywords internal
combine_chromosome_variant_data <- function(chr_results) {
  # Initialize combined data structure
  combined <- list(
    caqtl_data = list(),
    eqtl_data = list(),
    caqtl_indexed = list(),
    eqtl_indexed = list(),
    peak_gene_links = NULL,
    variant_peak_overlaps = list(),
    variant_metadata = list(),
    all_variants = character(0),
    n_variants = 0
  )

  # Get cell types from first chromosome result
  first_result <- chr_results[[1]]
  if (is.null(first_result)) {
    return(combined)
  }

  # Handle both naming conventions
  if (!is.null(first_result$caqtl_data_by_ct)) {
    cell_types <- names(first_result$caqtl_data_by_ct)
  } else if (!is.null(first_result$caqtl_data)) {
    cell_types <- names(first_result$caqtl_data)
  } else {
    cell_types <- character(0)
  }

  # Initialize collectors for rbindlist (avoids repeated rbind reallocation)
  caqtl_chunks <- list()
  eqtl_chunks <- list()
  for (ct in cell_types) {
    caqtl_chunks[[ct]] <- list()
    eqtl_chunks[[ct]] <- list()
    combined$caqtl_indexed[[ct]] <- list()
    combined$eqtl_indexed[[ct]] <- list()
  }

  # Collect data from each chromosome
  for (chr_data in chr_results) {
    if (!is.null(chr_data)) {
      # Handle both naming conventions
      caqtl_data_field <- if (!is.null(chr_data$caqtl_data_by_ct)) "caqtl_data_by_ct" else "caqtl_data"
      eqtl_data_field <- if (!is.null(chr_data$eqtl_data_by_ct)) "eqtl_data_by_ct" else "eqtl_data"

      for (ct in names(chr_data[[caqtl_data_field]])) {
        if (nrow(chr_data[[caqtl_data_field]][[ct]]) > 0) {
          caqtl_chunks[[ct]][[length(caqtl_chunks[[ct]]) + 1L]] <- chr_data[[caqtl_data_field]][[ct]]
        }
        if (length(chr_data$caqtl_indexed[[ct]]) > 0) {
          combined$caqtl_indexed[[ct]] <- c(combined$caqtl_indexed[[ct]], chr_data$caqtl_indexed[[ct]])
        }
      }

      for (ct in names(chr_data[[eqtl_data_field]])) {
        if (nrow(chr_data[[eqtl_data_field]][[ct]]) > 0) {
          eqtl_chunks[[ct]][[length(eqtl_chunks[[ct]]) + 1L]] <- chr_data[[eqtl_data_field]][[ct]]
        }
        if (length(chr_data$eqtl_indexed[[ct]]) > 0) {
          combined$eqtl_indexed[[ct]] <- c(combined$eqtl_indexed[[ct]], chr_data$eqtl_indexed[[ct]])
        }
      }

      combined$all_variants <- unique(c(combined$all_variants, chr_data$all_variants))
      combined$variant_peak_overlaps <- c(combined$variant_peak_overlaps, chr_data$variant_peak_overlaps)
      combined$variant_metadata <- c(combined$variant_metadata, chr_data$variant_metadata)
    }

    gc(verbose = FALSE)
  }

  # Combine all chunks at once using rbindlist (single allocation per cell type)
  for (ct in cell_types) {
    combined$caqtl_data[[ct]] <- if (length(caqtl_chunks[[ct]]) > 0) {
      data.table::rbindlist(caqtl_chunks[[ct]], fill = TRUE)
    } else {
      data.table::data.table()
    }
    combined$eqtl_data[[ct]] <- if (length(eqtl_chunks[[ct]]) > 0) {
      data.table::rbindlist(eqtl_chunks[[ct]], fill = TRUE)
    } else {
      data.table::data.table()
    }
  }

  # Use the peak-gene links from any chromosome (they should be the same)
  combined$peak_gene_links <- first_result$peak_gene_links

  # Update variant count
  combined$n_variants <- length(combined$all_variants)

  # Final garbage collection
  gc(verbose = FALSE)

  # Return with consistent field names that match load_variant_data_by_qtl_type
  return(list(
    all_variants = combined$all_variants,
    caqtl_data_by_ct = combined$caqtl_data,
    eqtl_data_by_ct = combined$eqtl_data,
    caqtl_indexed = combined$caqtl_indexed,
    eqtl_indexed = combined$eqtl_indexed,
    peak_gene_links = combined$peak_gene_links,
    variant_peak_overlaps = combined$variant_peak_overlaps,
    peak_data = combined$peak_data,
    n_variants = combined$n_variants
  ))
}

#' Load Variant Data (Two-Level)
#'
#' Load variant data with QTL information for two-level approach
#'
#' @param caqtl_files List of caQTL files by cell type
#' @param eqtl_files List of eQTL files by cell type
#' @param peak_gene_files Named list mapping cell types to peak-gene link file paths
#' @param chromosomes Chromosomes to analyze
#' @param pip_threshold Maximum PIP threshold across cell types for filtering (default: 0.5)
#' @param min_pip_threshold Minimum PIP threshold per cell type (default: 0.1)
#' @param acat_fdr_threshold FDR threshold for ACAT filtering
#' @param peak_bed_file Path to peak BED file
#' @param column_mapping Column name mapping
#' @param config Configuration list (used for cache and column mapping lookup)
#' @param quiet If TRUE, suppress informational logging
#' @param preloaded_peak_gene_data Optional pre-loaded peak-gene link tables, keyed by cell type
#' @return List with variant data organized by cell type
#' @noRd
load_variant_data_by_qtl_type <- function(caqtl_files, eqtl_files, peak_gene_files,
                                          chromosomes = NULL, pip_threshold = 0.5,
                                          min_pip_threshold = 0.1,
                                          acat_fdr_threshold = 0.05, peak_bed_file = NULL,
                                          column_mapping = NULL, config = NULL,
                                          quiet = FALSE,
                                          preloaded_peak_gene_data = NULL) {
  # Suppress logging if in parallel mode (quiet = TRUE)
  if (!quiet) {
    cli::cli_alert_info("Loading variant data for two-level analysis")
    cli::cli_alert_info("Hybrid PIP filtering:")
    cli::cli_alert_info("  - Max PIP threshold across cell types: {.val {pip_threshold}} (exclusive)")
    cli::cli_alert_info("  - Min PIP threshold per cell type: {.val {min_pip_threshold}} (exclusive)")
  }

  all_variants <- character(0)

  # Get column names from mapping or use defaults
  if (is.null(column_mapping)) {
    caqtl_variant_col <- "rsid"
    caqtl_pip_col <- "prob"
    caqtl_chr_col <- "chromosome"
    caqtl_trait_col <- "region"
    eqtl_variant_col <- "rsid"
    eqtl_pip_col <- "prob"
    eqtl_chr_col <- "chromosome"
    eqtl_trait_col <- "region"
    peak_gene_peak_col <- "peak_id"
    peak_gene_gene_col <- "gene_id"
  } else {
    caqtl_variant_col <- column_mapping$caqtl_susie$variant_id
    caqtl_pip_col <- column_mapping$caqtl_susie$pip
    caqtl_chr_col <- column_mapping$caqtl_susie$chromosome
    caqtl_trait_col <- column_mapping$caqtl_susie$feature_id
    eqtl_variant_col <- column_mapping$eqtl_susie$variant_id
    eqtl_pip_col <- column_mapping$eqtl_susie$pip
    eqtl_chr_col <- column_mapping$eqtl_susie$chromosome
    eqtl_trait_col <- column_mapping$eqtl_susie$feature_id
    peak_gene_peak_col <- column_mapping$peak_gene_links$peak_id %||% "peak_id"
    peak_gene_gene_col <- column_mapping$peak_gene_links$gene_id %||% "gene_id"
  }

  # Phase 1: Load all caQTL data with minimum PIP threshold to identify variants meeting max PIP criterion
  if (!quiet) cli::cli_alert_info("Phase 1: Loading caQTL data and calculating max PIPs")

  # Only read the essential columns (4 of 66) to reduce parse/memory by ~94%
  caqtl_select_cols <- c(caqtl_variant_col, caqtl_pip_col, caqtl_chr_col, caqtl_trait_col)

  # First pass: load all data with minimum threshold
  all_caqtl_data <- list()
  for (ct in names(caqtl_files)) {
    if (file.exists(caqtl_files[[ct]])) {
      if (!quiet) cli::cli_alert_info("Loading {ct} from {.file {caqtl_files[[ct]]}}")

      data <- data.table::data.table()
      tryCatch(
        {
          # Read only essential columns from SuSiE file
          data <- data.table::fread(caqtl_files[[ct]], select = caqtl_select_cols)

          # Standardize column names immediately after loading
          caqtl_rename <- list(
            variant_id = caqtl_variant_col, feature_id = caqtl_trait_col,
            pip = caqtl_pip_col, chromosome = caqtl_chr_col
          )
          rename_columns(data, caqtl_rename)

          # Apply minimum threshold only
          data <- data[pip > min_pip_threshold, ]
        },
        error = function(e) {
          if (!quiet) cli::cli_warn("Error loading {ct}: {e$message}")
        }
      )
      all_caqtl_data[[ct]] <- data
    } else {
      if (!quiet) cli::cli_warn("caQTL file not found: {.file {caqtl_files[[ct]]}}")
      all_caqtl_data[[ct]] <- data.table::data.table()
    }
  }

  # Calculate max PIP per variant-feature pair across all cell types
  all_caqtl_pips <- data.table::rbindlist(
    lapply(names(all_caqtl_data), function(ct) {
      data <- all_caqtl_data[[ct]]
      if (nrow(data) > 0 && "pip" %in% names(data)) {
        data.table::data.table(
          variant_id = data[["variant_id"]],
          feature_id = data[["feature_id"]],
          pip = data[["pip"]],
          cell_type = ct
        )
      }
    }),
    fill = TRUE
  )

  # Identify variant-feature pairs meeting max PIP threshold
  caqtl_keep_pairs <- data.table::data.table()
  if (nrow(all_caqtl_pips) > 0) {
    # Group by variant-feature pair to get max PIP across cell types
    max_caqtl_pips <- all_caqtl_pips[, .(max_pip = max(pip, na.rm = TRUE)), by = .(variant_id, feature_id)]
    caqtl_keep_pairs <- max_caqtl_pips[max_pip > pip_threshold, .(variant_id, feature_id)]

    # Create unique list of variants that have at least one qualifying feature
    caqtl_keep_variants <- unique(caqtl_keep_pairs$variant_id)

    if (!quiet) {
      cli::cli_alert_success("Identified {.val {nrow(caqtl_keep_pairs)}} caQTL variant-peak pairs with max PIP > {.val {pip_threshold}}")
      cli::cli_alert_success("  Involving {.val {length(caqtl_keep_variants)}} unique variants")
    }
  } else {
    caqtl_keep_variants <- character(0)
  }

  # Add qualifying caQTL variants to all_variants
  all_variants <- unique(c(all_variants, caqtl_keep_variants))

  # Phase 2: Filter each cell type's data to keep only qualified variant-feature pairs
  if (!quiet) cli::cli_alert_info("Phase 2: Applying hybrid filtering to caQTL data")
  caqtl_data_by_ct <- list()
  caqtl_indexed <- list()

  for (ct in names(all_caqtl_data)) {
    data <- all_caqtl_data[[ct]]
    if (nrow(data) > 0) {
      # Keep only variant-feature pairs meeting max PIP criterion (keyed join)
      data <- data[caqtl_keep_pairs, on = c("variant_id", "feature_id"), nomatch = NULL]

      if (!quiet && nrow(data) > 0) {
        cli::cli_alert_info("{ct}: {.val {nrow(data)}} variant-peak pairs retained after hybrid filtering")
      }

      caqtl_data_by_ct[[ct]] <- data

      # Create indexed version for lookup (only if data is not too large)
      if (nrow(data) < 50000) {
        caqtl_indexed[[ct]] <- split(data, data[["variant_id"]])
      } else {
        caqtl_indexed[[ct]] <- list()
        if (!quiet) cli::cli_alert_info("Skipping indexing for {ct} (too large: {nrow(data)} rows)")
      }
    } else {
      caqtl_data_by_ct[[ct]] <- data.table::data.table()
    }
  }

  # Phase 3: Load all eQTL data with minimum PIP threshold
  if (!quiet) cli::cli_alert_info("Phase 3: Loading eQTL data and calculating max PIPs")

  # Only read the essential columns (4 of 66) to reduce parse/memory by ~94%
  eqtl_select_cols <- c(eqtl_variant_col, eqtl_pip_col, eqtl_chr_col, eqtl_trait_col)

  # First pass: load all data with minimum threshold
  all_eqtl_data <- list()
  for (ct in names(eqtl_files)) {
    if (file.exists(eqtl_files[[ct]])) {
      if (!quiet) cli::cli_alert_info("Loading {ct} from {.file {eqtl_files[[ct]]}}")

      data <- data.table::data.table()
      tryCatch(
        {
          # Read only essential columns from SuSiE file
          data <- data.table::fread(eqtl_files[[ct]], select = eqtl_select_cols)

          # Standardize column names immediately after loading
          eqtl_rename <- list(
            variant_id = eqtl_variant_col, feature_id = eqtl_trait_col,
            pip = eqtl_pip_col, chromosome = eqtl_chr_col
          )
          rename_columns(data, eqtl_rename)

          # Apply minimum threshold only
          data <- data[pip > min_pip_threshold, ]
        },
        error = function(e) {
          if (!quiet) cli::cli_warn("Error loading {ct}: {e$message}")
        }
      )
      all_eqtl_data[[ct]] <- data
    } else {
      if (!quiet) cli::cli_warn("eQTL file not found: {.file {eqtl_files[[ct]]}}")
      all_eqtl_data[[ct]] <- data.table::data.table()
    }
  }

  # Calculate max PIP per variant-feature pair across all cell types
  all_eqtl_pips <- data.table::rbindlist(
    lapply(names(all_eqtl_data), function(ct) {
      data <- all_eqtl_data[[ct]]
      if (nrow(data) > 0 && "pip" %in% names(data)) {
        data.table::data.table(
          variant_id = data[["variant_id"]],
          feature_id = data[["feature_id"]],
          pip = data[["pip"]],
          cell_type = ct
        )
      }
    }),
    fill = TRUE
  )

  # Identify variant-feature pairs meeting max PIP threshold
  eqtl_keep_pairs <- data.table::data.table()
  if (nrow(all_eqtl_pips) > 0) {
    # Group by variant-feature pair to get max PIP across cell types
    max_eqtl_pips <- all_eqtl_pips[, .(max_pip = max(pip, na.rm = TRUE)), by = .(variant_id, feature_id)]
    eqtl_keep_pairs <- max_eqtl_pips[max_pip > pip_threshold, .(variant_id, feature_id)]

    # Create unique list of variants that have at least one qualifying feature
    eqtl_keep_variants <- unique(eqtl_keep_pairs$variant_id)

    if (!quiet) {
      cli::cli_alert_success("Identified {.val {nrow(eqtl_keep_pairs)}} eQTL variant-gene pairs with max PIP > {.val {pip_threshold}}")
      cli::cli_alert_success("  Involving {.val {length(eqtl_keep_variants)}} unique variants")
    }
  } else {
    eqtl_keep_variants <- character(0)
  }

  # Add qualifying eQTL variants to all_variants
  all_variants <- unique(c(all_variants, eqtl_keep_variants))

  # Phase 4: Filter each cell type's data to keep only qualified variant-feature pairs
  if (!quiet) cli::cli_alert_info("Phase 4: Applying hybrid filtering to eQTL data")
  eqtl_data_by_ct <- list()
  eqtl_indexed <- list()

  for (ct in names(all_eqtl_data)) {
    data <- all_eqtl_data[[ct]]
    if (nrow(data) > 0) {
      # Keep only variant-feature pairs meeting max PIP criterion (keyed join)
      data <- data[eqtl_keep_pairs, on = c("variant_id", "feature_id"), nomatch = NULL]

      if (!quiet && nrow(data) > 0) {
        cli::cli_alert_info("{ct}: {.val {nrow(data)}} variant-gene pairs retained after hybrid filtering")
      }

      eqtl_data_by_ct[[ct]] <- data

      # Create indexed version for lookup (only if data is not too large)
      if (nrow(data) < 50000) {
        eqtl_indexed[[ct]] <- split(data, data[["variant_id"]])
      } else {
        eqtl_indexed[[ct]] <- list()
        if (!quiet) cli::cli_alert_info("Skipping indexing for {ct} (too large: {nrow(data)} rows)")
      }
    } else {
      eqtl_data_by_ct[[ct]] <- data.table::data.table()
    }
  }

  # Load peak-gene link data
  if (!quiet) cli::cli_alert_info("Loading peak-gene link data")
  peak_gene_links <- NULL

  if (!is.null(peak_gene_files)) {
    if (!is.list(peak_gene_files) || is.null(names(peak_gene_files))) {
      cli::cli_abort("{.arg peak_gene_files} must be a named list mapping cell types to file paths.")
    }

    pg_list <- list()
    for (ct in names(peak_gene_files)) {
      # Use pre-loaded data if available (avoids re-reading same files per chromosome)
      if (!is.null(preloaded_peak_gene_data) && ct %in% names(preloaded_peak_gene_data)) {
        pg_data <- data.table::copy(preloaded_peak_gene_data[[ct]])

        # Rename columns to standard names if needed
        if (peak_gene_peak_col != "peak_id" && peak_gene_peak_col %in% names(pg_data)) {
          data.table::setnames(pg_data, peak_gene_peak_col, "peak_id")
        }
        if (peak_gene_gene_col != "gene_id" && peak_gene_gene_col %in% names(pg_data)) {
          data.table::setnames(pg_data, peak_gene_gene_col, "gene_id")
        }

        pg_data$config_cell_type <- ct
        pg_list[[ct]] <- pg_data
      } else if (file.exists(peak_gene_files[[ct]])) {
        pg_data <- data.table::fread(peak_gene_files[[ct]])

        # Validate required columns exist before renaming
        if (!all(c(peak_gene_peak_col, peak_gene_gene_col) %in% names(pg_data))) {
          cli::cli_abort("Peak-gene link file must contain columns '{peak_gene_peak_col}' and '{peak_gene_gene_col}'. Found: {paste(names(pg_data), collapse=', ')} in file: {peak_gene_files[[ct]]}")
        }

        # Rename columns to standard names if needed
        if (peak_gene_peak_col != "peak_id" || peak_gene_gene_col != "gene_id") {
          data.table::setnames(pg_data,
            old = c(peak_gene_peak_col, peak_gene_gene_col),
            new = c("peak_id", "gene_id")
          )
        }

        n_links <- nrow(pg_data)
        if (!quiet) {
          cli::cli_alert_info("Loaded {n_links} pre-filtered peak-gene links for {ct}")
        }

        pg_data$config_cell_type <- ct
        pg_list[[ct]] <- pg_data
      } else {
        cli::cli_abort("Peak-gene link file not found: {.file {peak_gene_files[[ct]]}}")
      }
    }

    # Combine all loaded data
    if (length(pg_list) > 0) {
      peak_gene_links <- data.table::rbindlist(pg_list, fill = TRUE)
    }
  }

  # Load peak data for overlap detection
  peak_data <- NULL
  variant_peak_overlaps <- list()

  if (is.null(peak_bed_file)) {
    cli::cli_abort("Peak BED file must be specified. Peak data is required for variant-peak overlap detection.")
  }

  # Handle cell-type specific peak_bed files
  if (grepl("\\{CELL_TYPE\\}", peak_bed_file)) {
    # Specification 1: Cell-type specific peak files
    if (!quiet) cli::cli_alert_info("Loading cell-type specific peak files")
    peak_data <- list()

    # Load peaks for each cell type
    for (ct in names(caqtl_files)) {
      ct_peak_file <- gsub("\\{CELL_TYPE\\}", ct, peak_bed_file)
      if (file.exists(ct_peak_file)) {
        if (!quiet) cli::cli_alert_info("Loading peaks for {ct}")
        peak_data[[ct]] <- load_peaks(ct_peak_file, target_chrs = chromosomes, quiet = quiet)
      } else {
        cli::cli_abort("Peak file not found for {ct}: {.file {ct_peak_file}}")
      }
    }
  } else if (file.exists(peak_bed_file)) {
    # Specification 2: Single peak file for all cell types
    if (!quiet) cli::cli_alert_info("Loading shared peak file for all cell types")
    peak_data <- load_peaks(peak_bed_file, target_chrs = chromosomes, quiet = quiet)
  } else {
    cli::cli_abort("Peak file not found: {.file {peak_bed_file}}")
  }

  # Calculate variant-peak overlaps
  if (!is.null(peak_data)) {
    if (!quiet) cli::cli_alert_info("Calculating variant-peak overlaps")

    # Debug: Check peak_data structure
    if (is.list(peak_data)) {
      if ("gr" %in% names(peak_data)) {
        if (!quiet) cli::cli_alert_info("Peak data structure: Single peak file with GRanges")
      } else if (length(peak_data) > 0 && !is.null(names(peak_data))) {
        if (!quiet) {
          cli::cli_alert_info("Peak data structure: Cell-type specific peaks")
          cli::cli_alert_info("Cell types found: {.val {names(peak_data)}}")
        }
      } else {
        if (!quiet) cli::cli_alert_info("Peak data structure: Unknown list format")
      }
    }

    # Get unique variant positions using vectorized operations
    if (length(all_variants) > 0) {
      # Use data.table for parsing
      var_dt <- data.table::data.table(variant_id = all_variants)
      # Split variant IDs vectorized
      var_parts <- data.table::tstrsplit(var_dt$variant_id, "_", fixed = TRUE)

      # Only keep variants with at least chr and pos
      valid_vars <- !is.na(var_parts[[1]]) & !is.na(var_parts[[2]])

      var_positions <- data.frame(
        variant_id = var_dt$variant_id[valid_vars],
        chr = var_parts[[1]][valid_vars],
        pos = as.numeric(var_parts[[2]][valid_vars]),
        stringsAsFactors = FALSE
      )
    } else {
      var_positions <- NULL
    }

    if (!is.null(var_positions) && nrow(var_positions) > 0) {
      # Create GRanges for variants
      gr_variants <- GenomicRanges::GRanges(
        seqnames = var_positions$chr,
        ranges = IRanges::IRanges(start = var_positions$pos, width = 1),
        variant_id = var_positions$variant_id
      )

      # Handle different peak data structures
      # Check if it's a single peak data structure (has $gr at top level)
      if (is.list(peak_data) && "gr" %in% names(peak_data) && inherits(peak_data$gr, "GRanges")) {
        # Single peak file for all cell types (specification 2)
        if (!quiet) cli::cli_alert_info("Using single peak file for overlap detection")
        tryCatch(
          {
            overlaps <- GenomicRanges::findOverlaps(gr_variants, peak_data$gr)
          },
          error = function(e) {
            cli::cli_alert_danger("ERROR in findOverlaps (single peak file):")
            cli::cli_alert_danger("gr_variants class: {.cls {class(gr_variants)}}")
            cli::cli_alert_danger("peak_data$gr class: {.cls {class(peak_data$gr)}}")
            cli::cli_alert_danger("Error message: {e$message}")
            stop(e)
          }
        )

        # Store overlaps by variant - for single peak file, no cell type info
        if (length(overlaps) > 0) {
          var_indices <- S4Vectors::queryHits(overlaps)
          peak_indices <- S4Vectors::subjectHits(overlaps)
          var_ids <- gr_variants$variant_id[var_indices]
          peak_ids <- peak_data$gr$peak_id[peak_indices]

          # For single peak file, just store peak IDs as before
          overlap_df <- data.frame(
            var_id = var_ids,
            peak_id = peak_ids,
            stringsAsFactors = FALSE
          )

          # Split by variant ID to create the list
          variant_peak_overlaps <- split(overlap_df$peak_id, overlap_df$var_id)
        }
      } else if (is.list(peak_data) && length(peak_data) > 0 && !is.null(names(peak_data))) {
        # Cell-type specific peak files (specification 1)
        # Check if this is indeed cell-type specific data
        first_ct <- names(peak_data)[1]
        if (!is.null(first_ct) && is.list(peak_data[[first_ct]]) && "gr" %in% names(peak_data[[first_ct]])) {
          if (!quiet) {
            cli::cli_alert_info("Using cell-type specific peak files for overlap detection")
            # Process peaks for each cell type separately to maintain cell type specificity
            cli::cli_alert_info("Processing cell-type specific peak overlaps")
          }

          for (ct in names(peak_data)) {
            if (!is.null(peak_data[[ct]]) &&
              is.list(peak_data[[ct]]) &&
              "gr" %in% names(peak_data[[ct]]) &&
              inherits(peak_data[[ct]]$gr, "GRanges")) {
              if (!quiet) cli::cli_alert_info("Processing {ct} peaks")
              ct_gr <- peak_data[[ct]]$gr

              # Find overlaps for this cell type
              overlaps <- GenomicRanges::findOverlaps(gr_variants, ct_gr)

              if (length(overlaps) > 0) {
                var_indices <- S4Vectors::queryHits(overlaps)
                peak_indices <- S4Vectors::subjectHits(overlaps)
                var_ids <- gr_variants$variant_id[var_indices]
                peak_ids <- ct_gr$peak_id[peak_indices]

                # Build overlap entries in one pass and group by variant
                overlap_entries <- Map(function(p) list(cell_type = ct, peak_id = p), peak_ids)
                names(overlap_entries) <- NULL
                overlap_by_var <- split(overlap_entries, var_ids)
                for (var_id in names(overlap_by_var)) {
                  variant_peak_overlaps[[var_id]] <- c(
                    variant_peak_overlaps[[var_id]],
                    overlap_by_var[[var_id]]
                  )
                }
              }
            }
          }
        } else {
          if (!quiet) cli::cli_warn("Peak data structure not recognized for overlap detection")
        }
      } else {
        if (!quiet) cli::cli_warn("Peak data is not in expected format for variant-peak overlap detection")
      }
    }
  }

  # Remove duplicates
  all_variants <- unique(all_variants)

  # Apply debug limit if enabled
  if (!is.null(config) && !is.null(config$debug) && config$debug$enabled) {
    # Use the same limit as peaks since variants are tied to peaks
    n_variants_limit <- config$debug$n_peaks
    if (length(all_variants) > n_variants_limit) {
      if (!quiet) cli::cli_alert_warning("Limiting to {.val {n_variants_limit}} variants out of {.val {length(all_variants)}}")
      all_variants <- all_variants[1:n_variants_limit]
    }
  }

  # Report hybrid filtering summary
  if (!quiet) {
    cli::cli_alert_success("Hybrid filtering summary:")
    cli::cli_alert_success("  - Total unique variants: {.val {length(all_variants)}}")
    cli::cli_alert_success("  - Retained variant-feature pairs with max PIP > {.val {pip_threshold}} across cell types")
    cli::cli_alert_success("  - And cell-type-specific PIP > {.val {min_pip_threshold}}")
  }

  return(list(
    all_variants = all_variants,
    caqtl_data_by_ct = caqtl_data_by_ct,
    eqtl_data_by_ct = eqtl_data_by_ct,
    caqtl_indexed = caqtl_indexed,
    eqtl_indexed = eqtl_indexed,
    peak_gene_links = peak_gene_links,
    variant_peak_overlaps = variant_peak_overlaps,
    peak_data = peak_data
  ))
}

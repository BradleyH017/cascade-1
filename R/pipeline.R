#' Run the CASCADE pipeline
#'
#' Main entry point for a CASCADE analysis. Loads configured QTL inputs,
#' runs gene/peak/variant categorization, and writes results to disk.
#'
#' @param config Configuration object (from `create_config()`) or path to a JSON config file.
#' @param output_dir Directory to write results to. Created if it does not exist.
#' @return Invisibly, a list of categorization results (gene, peak, variant).
#' @import cli
#' @import data.table
#' @export
run_cascade <- function(config, output_dir = "results") {
  # Load config if path provided
  if (is.character(config)) {
    cli::cli_alert_info("Loading configuration from: {.file {config}}")
    config <- jsonlite::fromJSON(config, simplifyVector = FALSE)
  }

  # Ensure config is a list
  if (!is.list(config)) {
    stop("Config must be a list object after loading, not ", class(config))
  }

  # Resolve cell type hierarchy (from JSON block, CellTypeHierarchy object, or default)
  config <- resolve_hierarchy(config)

  # Ensure column_mapping has defaults merged in
  if (is.null(config$column_mapping)) {
    config$column_mapping <- DEFAULT_COLUMN_MAPPING
  } else {
    for (file_type in names(DEFAULT_COLUMN_MAPPING)) {
      if (is.null(config$column_mapping[[file_type]])) {
        config$column_mapping[[file_type]] <- DEFAULT_COLUMN_MAPPING[[file_type]]
      }
    }
  }

  # Log Cochran's Q threshold from config if specified
  if (!is.null(config$parameters$cochran_q_threshold)) {
    cli::cli_alert_info("Cochran's Q threshold: {config$parameters$cochran_q_threshold}")
  } else {
    # Set default if not specified (genome-wide significance level)
    config$parameters$cochran_q_threshold <- 5e-8
    cli::cli_alert_info("Using default Cochran's Q threshold: 5e-8")
  }

  # Set LFSR thresholds from config or use defaults
  if (is.null(config$parameters$lfsr_sig_threshold)) {
    config$parameters$lfsr_sig_threshold <- LFSR_SIG_THRESHOLD
    cli::cli_alert_info("Using default LFSR significance threshold: {LFSR_SIG_THRESHOLD}")
  } else {
    cli::cli_alert_info("LFSR significance threshold: {config$parameters$lfsr_sig_threshold}")
  }

  if (is.null(config$parameters$lfsr_null_threshold)) {
    config$parameters$lfsr_null_threshold <- LFSR_NULL_THRESHOLD
    cli::cli_alert_info("Using default LFSR null threshold: {LFSR_NULL_THRESHOLD}")
  } else {
    cli::cli_alert_info("LFSR null threshold: {config$parameters$lfsr_null_threshold}")
  }

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Always analyze all feature types
  feature_types <- c("gene", "peak", "variant")
  cli::cli_alert_info("Running comprehensive analysis for all feature types: gene, peak, and variant")

  results <- list()

  # Step 1: Load data for all analyses
  cli::cli_h2("Loading Data")

  # Load gene data
  cli::cli_alert("Loading gene data...")
  gene_data <- load_feature_data(config, "gene", config$chromosomes,
    num_cores = config$parameters$n_cores
  )

  # Load peak data
  cli::cli_alert("Loading peak data...")
  peak_data <- load_feature_data(config, "peak", config$chromosomes,
    num_cores = config$parameters$n_cores
  )

  # Load variant data
  cli::cli_alert("Loading variant data...")
  variant_data <- load_variant_data(
    config = config,
    chromosomes = config$chromosomes,
    pip_threshold = config$parameters$pip_threshold,
    min_pip_threshold = config$parameters$min_pip_threshold %||% 0.1,
    acat_fdr_threshold = config$parameters$acat_fdr_threshold,
    peak_bed_file = config$file_patterns$peak_bed,
    column_mapping = config$column_mapping,
    num_cores = config$parameters$n_cores
  )

  # Step 2: Load LFSR results and mashr models if available
  cli::cli_h2("Loading Pre-computed LFSR Results")
  lfsr_results <- load_lfsr_results(config)

  # Load pre-computed meta data with Cochran's Q values
  cli::cli_h2("Loading Pre-computed Meta Data")
  meta_data <- list()

  # Initialize cache once for meta data loading
  cache_obj <- NULL
  if (!is.null(config$cache) && !is.null(config$cache$enabled) && config$cache$enabled) {
    cache_obj <- init_cache(config)
  }

  # Load required meta data for both gene and peak analyses
  # eQTL meta data is required
  if (is.null(config$file_patterns$eqtl_meta)) {
    cli::cli_abort("eQTL meta data file is required but not specified in config. Please add 'eqtl_meta' to file_patterns in your configuration.")
  }
  validate_file_exists(config$file_patterns$eqtl_meta, "eQTL meta data file")
  cli::cli_alert_info("Loading eQTL meta data from: {.file {config$file_patterns$eqtl_meta}}")
  meta_data$eqtl <- load_meta_data(config$file_patterns$eqtl_meta, cache = cache_obj, column_mapping = config$column_mapping$meta)
  if (is.null(meta_data$eqtl) || nrow(meta_data$eqtl) == 0) {
    cli::cli_abort("Failed to load eQTL meta data or file is empty: {.file {config$file_patterns$eqtl_meta}}")
  }

  # caQTL meta data is required
  if (is.null(config$file_patterns$caqtl_meta)) {
    cli::cli_abort("caQTL meta data file is required but not specified in config. Please add 'caqtl_meta' to file_patterns in your configuration.")
  }
  validate_file_exists(config$file_patterns$caqtl_meta, "caQTL meta data file")
  cli::cli_alert_info("Loading caQTL meta data from: {.file {config$file_patterns$caqtl_meta}}")
  meta_data$caqtl <- load_meta_data(config$file_patterns$caqtl_meta, cache = cache_obj, column_mapping = config$column_mapping$meta)
  if (is.null(meta_data$caqtl) || nrow(meta_data$caqtl) == 0) {
    cli::cli_abort("Failed to load caQTL meta data or file is empty: {.file {config$file_patterns$caqtl_meta}}")
  }

  # Load CS cluster data if available for improved heterogeneity analysis
  cs_clusters <- NULL
  cs_cluster_variants <- NULL

  if (!is.null(config$file_patterns$cs_clusters) && !is.null(config$file_patterns$cs_cluster_variants)) {
    cli::cli_h2("Loading CS Cluster Data for Improved Heterogeneity Analysis")

    # Load CS cluster mappings
    if (file.exists(config$file_patterns$cs_clusters)) {
      cs_clusters <- load_cs_clusters(config$file_patterns$cs_clusters, cache = cache_obj, column_mapping = config$column_mapping$cs_clusters)
      if (is.null(cs_clusters)) {
        cli::cli_alert_warning("Failed to load CS clusters, will use standard heterogeneity analysis")
      }
    } else {
      cli::cli_alert_info("CS clusters file not found: {.file {config$file_patterns$cs_clusters}}")
    }

    # Load CS cluster variants
    if (file.exists(config$file_patterns$cs_cluster_variants)) {
      cs_cluster_variants <- load_cs_cluster_variants(config$file_patterns$cs_cluster_variants, cache = cache_obj, column_mapping = config$column_mapping$cs_cluster_variants)
      if (is.null(cs_cluster_variants)) {
        cli::cli_alert_warning("Failed to load CS cluster variants, will use standard heterogeneity analysis")
      }
    } else {
      cli::cli_alert_info("CS cluster variants file not found: {.file {config$file_patterns$cs_cluster_variants}}")
    }

    # Only use CS clusters if both files loaded successfully
    if (!is.null(cs_clusters) && !is.null(cs_cluster_variants)) {
      cli::cli_alert_success("CS cluster data loaded successfully - will use cluster-based heterogeneity analysis")
    } else {
      cs_clusters <- NULL
      cs_cluster_variants <- NULL
      cli::cli_alert_info("Falling back to standard Cochran's Q-based heterogeneity analysis")
    }
  } else {
    cli::cli_alert_info("CS cluster files not configured - using standard heterogeneity analysis")
  }

  # Load mashr models if specified (use cached_readRDS for qs2 binary caching)
  mashr_models <- list()
  if (!is.null(config$file_patterns$eqtl_mashr)) {
    validate_file_exists(config$file_patterns$eqtl_mashr, "eQTL mashr model file")
    cli::cli_alert_info("Loading eQTL mashr model from: {.file {config$file_patterns$eqtl_mashr}}")
    mashr_models$eqtl <- cached_readRDS(config$file_patterns$eqtl_mashr, cache = cache_obj)
  }
  if (!is.null(config$file_patterns$caqtl_mashr)) {
    validate_file_exists(config$file_patterns$caqtl_mashr, "caQTL mashr model file")
    cli::cli_alert_info("Loading caQTL mashr model from: {.file {config$file_patterns$caqtl_mashr}}")
    mashr_models$caqtl <- cached_readRDS(config$file_patterns$caqtl_mashr, cache = cache_obj)
  }

  # Log LFSR loading details
  cli::cli_alert_info("LFSR rows loaded: eQTL long={.val {nrow(lfsr_results$eqtl_long)}}, caQTL long={.val {nrow(lfsr_results$caqtl_long)}}, eQTL feature={.val {nrow(lfsr_results$eqtl_feature)}}, caQTL feature={.val {nrow(lfsr_results$caqtl_feature)}}")

  # Step 3: Variant categorization (independent, run first)
  cli::cli_h2("Running Variant Categorization")

  # Categorize variants with pre-computed LFSR results (always two-stage)
  variant_categories <- categorize_variants(
    variant_data,
    lfsr_results = lfsr_results,
    config = config
  )

  # Extract variant-feature specificity for use in gene/peak categorization
  variant_feature_specificity <- variant_categories$variant_feature_specificity


  # Save two-stage results
  cli::cli_alert_info("Saving two-stage variant categorization results...")

  # Save per-cell-type results (Stage 1)
  save_variant_results_per_celltype(
    variant_categories$per_celltype,
    output_dir
  )

  # Log per-cell-type QTL detection counts (patterns 23-25 = no QTL)
  for (ct in names(variant_categories$per_celltype)) {
    ct_data <- variant_categories$per_celltype[[ct]]
    n_qtl <- sum(!(ct_data$qtl_pattern_number %in% c(23L, 24L, 25L)))
    cli::cli_alert_info("Per-cell-type QTL effects in {.val {ct}}: {.val {n_qtl}} variants")
  }

  # Save cross-cell-type results (Stage 2 - Combined L1+L2)
  if (!is.null(variant_categories$cross_celltype) && nrow(variant_categories$cross_celltype) > 0) {
    save_variant_results_cross_celltype(
      variant_categories$cross_celltype,
      output_dir
    )

    # Report the total number of cell types analyzed
    n_total_types <- length(config$cell_types)
    cli::cli_alert_info("Cross-cell-type aggregation completed ({n_total_types} cell types)")
  } else {
    cli::cli_alert_info("No cross-cell-type results available")
  }

  # Store both stages in results
  results$variant_categories <- variant_categories

  # Step 4: Gene categorization
  cli::cli_h2("Running Gene Categorization")
  cli::cli_alert_info("Processing {.val {gene_data$n_features}} genes")


  # Use categorization that extracts variants during analysis
  gene_result <- categorize_features(
    gene_data,
    lfsr_results = lfsr_results,
    susie_results = gene_data$susie_results,
    meta_data = meta_data$eqtl, # Required: pre-computed meta data
    variant_feature_specificity = variant_feature_specificity, # Required: variant-feature specificity
    feature_type = "gene",
    hierarchy = config$hierarchy,
    lfsr_sig_threshold = config$parameters$lfsr_sig_threshold,
    lfsr_null_threshold = config$parameters$lfsr_null_threshold,
    cochran_q_threshold = config$parameters$cochran_q_threshold,
    use_cs_clusters = !is.null(cs_clusters) && !is.null(cs_cluster_variants),
    cs_clusters = cs_clusters,
    cs_cluster_variants = cs_cluster_variants
  )

  # Extract results and variant details
  gene_categories <- gene_result$categories
  gene_top_variants <- gene_result$variant_details

  # Save results
  output_file <- file.path(output_dir, "gene_categorization.tsv.gz")
  data.table::fwrite(gene_categories, output_file,
    sep = "\t",
    quote = FALSE, compress = "gzip", na = "NA"
  )
  cli::cli_alert_success("Gene categorization results saved to: {.file {output_file}}")

  # Save top variants (already extracted during categorization)
  if (!is.null(gene_top_variants) && nrow(gene_top_variants) > 0) {
    # Determine filename based on whether CS clusters were used
    use_cs_clusters <- !is.null(cs_clusters) && !is.null(cs_cluster_variants)

    if (use_cs_clusters) {
      # CS cluster variant details have different structure
      top_var_file <- file.path(output_dir, "gene_categorization.cs_cluster_variants.tsv.gz")
      cli::cli_alert_success("Gene CS cluster top variants saved to: {.file {top_var_file}}")
      cli::cli_alert_info("Extracted {.val {nrow(gene_top_variants)}} CS cluster variant records for {.val {length(unique(gene_top_variants$gene_id))}} genes")

      # Also extract and save per-cell-type SuSiE results for CS cluster variants
      cli::cli_alert_info("Extracting per-cell-type SuSiE results for CS cluster variants...")
      gene_susie_details <- extract_cs_cluster_susie_details(
        features = unique(gene_categories$gene_id),
        cs_clusters = cs_clusters,
        cs_cluster_variants = cs_cluster_variants,
        susie_results = gene_data$susie_results,
        feature_type = "gene"
      )

      if (!is.null(gene_susie_details) && nrow(gene_susie_details) > 0) {
        # Add cluster patterns from the top variants analysis
        cluster_patterns <- unique(gene_top_variants[, .(cluster_id, cluster_pattern)])
        gene_susie_details <- merge(gene_susie_details, cluster_patterns,
          by = "cluster_id", all.x = TRUE
        )

        susie_file <- file.path(output_dir, "gene_categorization.cs_cluster_susie_details.tsv.gz")
        data.table::fwrite(gene_susie_details, susie_file,
          sep = "\t",
          quote = FALSE, compress = "gzip", na = "NA"
        )
        cli::cli_alert_success("Gene CS cluster SuSiE details saved to: {.file {susie_file}}")
        cli::cli_alert_info("Extracted {.val {nrow(gene_susie_details)}} SuSiE records across {.val {length(unique(gene_susie_details$cell_type))}} cell types")
      }
    } else {
      # Standard top variants output
      top_var_file <- file.path(output_dir, "gene_categorization.top_variants.tsv.gz")
      cli::cli_alert_success("Gene top variants saved to: {.file {top_var_file}}")
      cli::cli_alert_info("Extracted {.val {nrow(gene_top_variants)}} top variant records for {.val {length(unique(gene_top_variants$gene_id))}} genes")
    }

    # Flatten any list columns to character strings before writing
    list_cols <- names(gene_top_variants)[sapply(gene_top_variants, is.list)]
    if (length(list_cols) > 0) {
      for (lc in list_cols) {
        data.table::set(gene_top_variants,
          j = lc,
          value = vapply(gene_top_variants[[lc]], function(x) {
            if (is.null(x) || length(x) == 0) {
              NA_character_
            } else {
              paste(x, collapse = ",")
            }
          }, character(1))
        )
      }
    }

    data.table::fwrite(gene_top_variants, top_var_file,
      sep = "\t",
      quote = FALSE, compress = "gzip", na = "NA"
    )
  }

  # Generate and save summary
  gene_summary <- summarize_column(gene_categories, "cell_type_specificity", "Cell Type Specificity")
  summary_file <- file.path(output_dir, "gene_categorization.summary.tsv")
  save_categorization_summary(gene_summary, summary_file)

  results$gene_categories <- gene_categories

  # Step 5: Peak categorization
  cli::cli_h2("Running Peak Categorization")
  cli::cli_alert_info("Processing {.val {peak_data$n_features}} peaks")

  # Use categorization that extracts variants during analysis
  peak_result <- categorize_features(
    peak_data,
    lfsr_results = lfsr_results,
    susie_results = peak_data$susie_results,
    meta_data = meta_data$caqtl, # Required: pre-computed meta data
    variant_feature_specificity = variant_feature_specificity, # Required: variant-feature specificity
    feature_type = "peak",
    hierarchy = config$hierarchy,
    lfsr_sig_threshold = config$parameters$lfsr_sig_threshold,
    lfsr_null_threshold = config$parameters$lfsr_null_threshold,
    cochran_q_threshold = config$parameters$cochran_q_threshold,
    use_cs_clusters = !is.null(cs_clusters) && !is.null(cs_cluster_variants),
    cs_clusters = cs_clusters,
    cs_cluster_variants = cs_cluster_variants
  )

  # Extract results and variant details
  peak_categories <- peak_result$categories
  peak_top_variants <- peak_result$variant_details

  # Save results
  output_file <- file.path(output_dir, "peak_categorization.tsv.gz")
  data.table::fwrite(peak_categories, output_file,
    sep = "\t",
    quote = FALSE, compress = "gzip", na = "NA"
  )
  cli::cli_alert_success("Peak categorization results saved to: {.file {output_file}}")

  # Save top variants (already extracted during categorization)
  if (!is.null(peak_top_variants) && nrow(peak_top_variants) > 0) {
    # Determine filename based on whether CS clusters were used
    use_cs_clusters <- !is.null(cs_clusters) && !is.null(cs_cluster_variants)

    if (use_cs_clusters) {
      # CS cluster variant details have different structure
      top_var_file <- file.path(output_dir, "peak_categorization.cs_cluster_variants.tsv.gz")
      cli::cli_alert_success("Peak CS cluster top variants saved to: {.file {top_var_file}}")
      cli::cli_alert_info("Extracted {.val {nrow(peak_top_variants)}} CS cluster variant records for {.val {length(unique(peak_top_variants$peak_id))}} peaks")

      # Also extract and save per-cell-type SuSiE results for CS cluster variants
      cli::cli_alert_info("Extracting per-cell-type SuSiE results for CS cluster variants...")
      peak_susie_details <- extract_cs_cluster_susie_details(
        features = unique(peak_categories$peak_id),
        cs_clusters = cs_clusters,
        cs_cluster_variants = cs_cluster_variants,
        susie_results = peak_data$susie_results,
        feature_type = "peak"
      )

      if (!is.null(peak_susie_details) && nrow(peak_susie_details) > 0) {
        # Add cluster patterns from the top variants analysis
        cluster_patterns <- unique(peak_top_variants[, .(cluster_id, cluster_pattern)])
        peak_susie_details <- merge(peak_susie_details, cluster_patterns,
          by = "cluster_id", all.x = TRUE
        )

        susie_file <- file.path(output_dir, "peak_categorization.cs_cluster_susie_details.tsv.gz")
        data.table::fwrite(peak_susie_details, susie_file,
          sep = "\t",
          quote = FALSE, compress = "gzip", na = "NA"
        )
        cli::cli_alert_success("Peak CS cluster SuSiE details saved to: {.file {susie_file}}")
        cli::cli_alert_info("Extracted {.val {nrow(peak_susie_details)}} SuSiE records across {.val {length(unique(peak_susie_details$cell_type))}} cell types")
      }
    } else {
      # Standard top variants output
      top_var_file <- file.path(output_dir, "peak_categorization.top_variants.tsv.gz")
      cli::cli_alert_success("Peak top variants saved to: {.file {top_var_file}}")
      cli::cli_alert_info("Extracted {.val {nrow(peak_top_variants)}} top variant records for {.val {length(unique(peak_top_variants$peak_id))}} peaks")
    }

    # Flatten any list columns to character strings before writing
    list_cols <- names(peak_top_variants)[sapply(peak_top_variants, is.list)]
    if (length(list_cols) > 0) {
      for (lc in list_cols) {
        data.table::set(peak_top_variants,
          j = lc,
          value = vapply(peak_top_variants[[lc]], function(x) {
            if (is.null(x) || length(x) == 0) {
              NA_character_
            } else {
              paste(x, collapse = ",")
            }
          }, character(1))
        )
      }
    }

    data.table::fwrite(peak_top_variants, top_var_file,
      sep = "\t",
      quote = FALSE, compress = "gzip", na = "NA"
    )
  }

  # Generate and save summary
  peak_summary <- summarize_column(peak_categories, "cell_type_specificity", "Cell Type Specificity")
  summary_file <- file.path(output_dir, "peak_categorization.summary.tsv")
  save_categorization_summary(peak_summary, summary_file)

  results$peak_categories <- peak_categories

  cli::cli_h2("Analysis Complete")
  invisible(results)
}

#' Create Configuration Object
#'
#' Helper function to create a configuration object for cascade analysis
#'
#' @param cell_types Vector of cell type names
#' @param chromosomes Vector of chromosomes to analyze
#' @param file_patterns List with file pattern templates
#' @param parameters Analysis parameters including:
#'   \itemize{
#'     \item pip_threshold: Maximum PIP threshold across cell types for additional filtering of 95% CS variants (default: 0.5)
#'     \item min_pip_threshold: Minimum PIP threshold per cell type (default: 0.1)
#'     \item acat_fdr_threshold: FDR threshold for ACAT significance (default: 0.05)
#'     \item lfsr_sig_threshold: LFSR significance threshold (default: 0.05)
#'     \item lfsr_null_threshold: LFSR null hypothesis threshold (default: 0.5)
#'     \item run_mash: Whether to run mash analysis (default: FALSE)
#'     \item n_cores: Number of cores for parallelization (default: NULL for auto-detect)
#'     \item mash_params: Parameters for mash analysis
#'   }
#' @param column_mapping Column name mappings for input files
#' @param feature_type Type of features to analyze ("gene", "peak", "variant", or "all")
#' @return Configuration list
#' @export
create_config <- function(cell_types,
                          chromosomes = NULL,
                          file_patterns = list(),
                          parameters = list(
                            pip_threshold = 0.5,
                            min_pip_threshold = 0.1,
                            acat_fdr_threshold = 0.05,
                            lfsr_sig_threshold = LFSR_SIG_THRESHOLD,
                            lfsr_null_threshold = LFSR_NULL_THRESHOLD,
                            run_mash = FALSE,
                            n_cores = NULL,
                            mash_params = list(
                              max_variants_per_gene = 5,
                              alpha = 1,
                              strong_z_threshold = 2
                            )
                          ),
                          column_mapping = list(),
                          feature_type = "all") {
  # Merge user column_mapping overrides with defaults (per file type)
  merged_mapping <- DEFAULT_COLUMN_MAPPING
  for (file_type in names(column_mapping)) {
    if (file_type %in% names(merged_mapping)) {
      merged_mapping[[file_type]] <- modifyList(
        merged_mapping[[file_type]], column_mapping[[file_type]]
      )
    } else {
      merged_mapping[[file_type]] <- column_mapping[[file_type]]
    }
  }

  config <- list(
    cell_types = cell_types,
    chromosomes = chromosomes,
    file_patterns = file_patterns,
    parameters = parameters,
    column_mapping = merged_mapping,
    feature_type = feature_type
  )

  # Parse feature types
  if (grepl(",", feature_type)) {
    feature_types <- trimws(strsplit(feature_type, ",")[[1]])
  } else {
    feature_types <- feature_type
  }

  # Validate required file patterns
  required_patterns <- c()
  if ("gene" %in% feature_types) {
    required_patterns <- c(required_patterns, "eqtl_acat")
  }
  if ("peak" %in% feature_types) {
    required_patterns <- c(required_patterns, "caqtl_acat")
  }
  if ("variant" %in% feature_types) {
    required_patterns <- c(required_patterns, "caqtl_susie", "eqtl_susie")
  }

  missing_patterns <- setdiff(required_patterns, names(file_patterns))
  if (length(missing_patterns) > 0) {
    cli::cli_abort("Missing required file patterns: {.val {missing_patterns}}")
  }

  return(config)
}

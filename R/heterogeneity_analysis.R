#' @title Variant Heterogeneity Analysis
#' @description Functions for analyzing variant heterogeneity patterns using
#'   Cochran's Q test and CS cluster-based approaches

#' Add Variant Analysis to Categorization Results
#'
#' Add variant heterogeneity analysis using Cochran's Q test or CS clusters
#'
#' @param results Results from vectorized categorization
#' @param susie_results SuSiE results
#' @param data_type Either "gene" or "peak"
#' @param lfsr_results Optional LFSR results
#' @param meta_data Required. Pre-computed meta data with Cochran's Q values
#' @param variant_feature_specificity Required. Variant-feature specificity mapping
#' @param cochran_q_threshold Cochran's Q p-value threshold
#' @param use_cs_clusters Logical. Whether to use CS cluster-based analysis (default TRUE)
#' @param cs_clusters Optional. CS cluster mappings from load_cs_clusters()
#' @param cs_cluster_variants Optional. CS cluster variant mappings from load_cs_cluster_variants()
#' @return List with updated results and variant details
add_variant_analysis <- function(results, susie_results, data_type, lfsr_results,
                                 meta_data, variant_feature_specificity,
                                 cochran_q_threshold = 5e-8,
                                 use_cs_clusters = TRUE,
                                 cs_clusters = NULL,
                                 cs_cluster_variants = NULL) {
  # Check if results is empty
  if (is.null(results) || nrow(results) == 0) {
    return(list(results = results, variant_details = NULL))
  }

  # Get feature column name
  feature_col <- if (data_type == "gene") "gene_id" else "peak_id"

  # Use CS cluster-based analysis if requested and data is available
  if (use_cs_clusters && !is.null(cs_clusters) && !is.null(cs_cluster_variants)) {
    cli::cli_alert_info("Using CS cluster-based heterogeneity analysis")

    # Get features to analyze
    features <- results[[feature_col]]

    # Perform CS cluster analysis with variant details
    cluster_analysis <- analyze_cs_cluster_heterogeneity(
      features = features,
      cs_clusters = cs_clusters,
      cs_cluster_variants = cs_cluster_variants,
      meta_data = meta_data,
      cochran_q_threshold = cochran_q_threshold,
      feature_type = data_type,
      return_variant_details = TRUE
    )

    # Extract results and variant details
    cluster_results <- cluster_analysis$results
    cluster_variant_details <- cluster_analysis$variant_details

    # Merge cluster results with existing results
    # The cluster_results already has feature_id column, need to match it with our feature column
    data.table::setnames(cluster_results, "feature_id", feature_col)

    # Properly merge by feature ID to ensure correct alignment
    # First, store the original column order
    orig_cols <- names(results)

    # Remove columns that will be updated from results
    cols_to_update <- c(
      "variant_heterogeneity_pattern", "n_cs_clusters",
      "top_variant", "max_chisq", "variant_heterogeneity", "n_variants"
    )
    cols_to_remove <- intersect(cols_to_update, names(results))
    if (length(cols_to_remove) > 0) {
      results[, (cols_to_remove) := NULL]
    }

    # Merge with cluster results
    results <- merge(results, cluster_results, by = feature_col, all.x = TRUE, sort = FALSE)

    # Use n_cs_clusters as proxy for n_variants if needed
    if ("n_variants" %in% orig_cols && !"n_variants" %in% names(results)) {
      results[, n_variants := n_cs_clusters]
    }

    # Reorder columns to match original order plus new columns
    new_cols <- setdiff(names(results), orig_cols)
    data.table::setcolorder(results, c(orig_cols, new_cols))

    # For CS cluster analysis, return cluster-specific variant details
    cli::cli_alert_success("CS cluster analysis complete")

    return(list(
      results = results,
      variant_details = cluster_variant_details
    ))
  }

  # Otherwise, use the original implementation
  cli::cli_alert_info("Using Cochran's Q-based heterogeneity analysis")

  spec_col <- "cell_type_specificity"


  # Set keys for efficient lookup
  data.table::setkey(variant_feature_specificity, variant_id, feature_id, feature_type)

  n_features <- nrow(results)
  feature_col <- if (data_type == "gene") "gene" else "peak"

  # Initialize new columns
  results$variant_heterogeneity <- rep(NA, n_features)
  results$variant_heterogeneity_pattern <- rep(NA, n_features)
  results$n_variants <- rep(0, n_features)
  results$cs_cts <- rep("", n_features)
  results$n_cs_cts <- rep(0, n_features)

  # Get feature ID column name
  feature_col <- names(results)[1]

  # Process all features at once with extraction
  cli::cli_alert_info("Processing variant patterns and extracting variant details")

  # Prepare data for C++ function
  feature_ids <- results[[feature_col]]
  sig_cts_list <- lapply(seq_len(n_features), function(i) {
    sig_cts <- results$significant_cts[i]
    if (is.na(sig_cts) || sig_cts == "") {
      character(0)
    } else {
      strsplit(sig_cts, ",")[[1]]
    }
  })

  # Use pre-computed meta data (always required); must have direction and max_pip columns.
  meta_filtered <- meta_data

  # Create feature data.table
  feature_dt <- data.table::data.table(
    feature_id = feature_ids,
    idx = seq_len(n_features),
    n_sig_cts = vapply(sig_cts_list, length, integer(1))
  )

  # Only process features with 2+ significant cell types
  feature_dt[, process := n_sig_cts >= 2]

  # Find best variant per feature using data.table aggregation
  # This is MUCH faster than looping in C++
  if (nrow(meta_filtered) > 0) {
    best_variants <- meta_filtered[
      feature_id %in% feature_ids,
      .(
        top_variant = variant_id[which.max(max_chisq)],
        max_pip = max(max_pip, na.rm = TRUE),
        max_chisq = max(max_chisq, na.rm = TRUE),
        cochran_q_pval = cochran_q_pval[which.max(max_chisq)],
        direction = direction[which.max(max_chisq)]
      ),
      by = feature_id
    ]
  } else {
    # No high-PIP variants found
    best_variants <- data.table::data.table(
      feature_id = character(),
      top_variant = character(),
      max_pip = numeric(),
      max_chisq = numeric(),
      cochran_q_pval = numeric(),
      direction = character()
    )
  }

  # Determine patterns using vectorized operations (only if we have variants)
  if (nrow(best_variants) > 0) {
    best_variants[, pattern := ifelse(
      is.na(cochran_q_pval), NA_character_,
      ifelse(
        cochran_q_pval < cochran_q_threshold,
        ifelse(direction == "+-", "shared_opposite", "shared_heterogeneous"),
        "shared_consistent"
      )
    )]
  }

  # Count all variants per feature (not just high PIP ones)
  variant_counts <- meta_data[, .(n_variants = .N), by = feature_id]

  # Join results back to features using data.table merge
  if (nrow(best_variants) > 0) {
    data.table::setkey(best_variants, feature_id)
  }

  if (nrow(variant_counts) > 0) {
    data.table::setkey(variant_counts, feature_id)
  }

  data.table::setkey(feature_dt, feature_id)

  # Merge all results using left joins
  if (nrow(best_variants) > 0) {
    feature_dt <- merge(feature_dt, best_variants, by = "feature_id", all.x = TRUE)
  }

  if (nrow(variant_counts) > 0) {
    feature_dt <- merge(feature_dt, variant_counts, by = "feature_id", all.x = TRUE)
  }

  # Handle missing values for features without variants
  feature_dt[!process | is.na(pattern), pattern := NA_character_]
  feature_dt[!process | is.na(max_pip), max_pip := NA_real_]
  feature_dt[!process | is.na(top_variant), top_variant := NA_character_]
  feature_dt[is.na(n_variants), n_variants := 0L]

  # Features with no high-PIP variants get "distinct_variants"
  feature_dt[process == TRUE & is.na(pattern), pattern := "distinct_variants"]

  # Sort back to original order
  data.table::setkey(feature_dt, idx)

  # Update results
  results$variant_heterogeneity_pattern <- feature_dt$pattern
  results$max_pip <- feature_dt$max_pip
  results$max_chisq <- feature_dt$max_chisq
  results$top_variant <- feature_dt$top_variant
  results$n_variants <- feature_dt$n_variants

  # Set heterogeneity flag for heterogeneous patterns
  heterogeneous_patterns <- c("shared_opposite", "shared_heterogeneous")
  results$variant_heterogeneity <- results$variant_heterogeneity_pattern %in% heterogeneous_patterns

  # Add hierarchical variant patterns using data.table operations
  if (nrow(meta_filtered) > 0) {
    cli::cli_alert_info("Computing hierarchical variant heterogeneity patterns")

    # Create feature data.table with necessary info
    # Use the correct column name for specificity
    spec_values <- results[[spec_col]]

    features_dt <- data.table::data.table(
      feature_id = feature_ids,
      idx = seq_len(n_features),
      n_significant_cts = results$n_significant_cts,
      specificity_category = spec_values,
      significant_cts_str = vapply(sig_cts_list, paste, character(1), collapse = ",")
    )

    # Filter to features with 2+ significant cell types
    features_to_process <- features_dt[n_significant_cts >= 2]

    if (nrow(features_to_process) > 0) {
      # Filter variant_feature_specificity for this data type
      vfs_filtered <- variant_feature_specificity[feature_type == data_type]

      # Join meta_filtered with variant_feature_specificity to get all info at once
      meta_with_spec <- merge(
        meta_filtered[, .(feature_id, variant_id, direction, cochran_q_pval)],
        vfs_filtered[, .(feature_id, variant_id, variant_specificity, significant_cell_types)],
        by = c("feature_id", "variant_id"),
        all.x = TRUE
      )

      # Check for missing specs - this is expected when meta data and SuSiE data are from different analyses
      missing_specs <- meta_with_spec[is.na(variant_specificity)]
      if (nrow(missing_specs) > 0) {
        pct_missing <- round(100 * nrow(missing_specs) / nrow(meta_with_spec), 1)

        # Only show a note if some pairs matched (otherwise it's expected)
        if (pct_missing < 100) {
          cli::cli_alert_info("Note: {nrow(missing_specs)} ({pct_missing}%) variant-feature pairs from meta data don't have variant specificity from SuSiE")
        }

        # For missing specs, set a default value
        meta_with_spec[is.na(variant_specificity), variant_specificity := "Not available"]
        meta_with_spec[is.na(significant_cell_types), significant_cell_types := ""]
      }

      # Join with features to get feature-level info
      variant_analysis <- merge(
        meta_with_spec,
        features_to_process[, .(feature_id, specificity_category, significant_cts_str)],
        by = "feature_id"
      )

      # Pre-split cell type strings once (avoids redundant per-row strsplit)
      var_cts_split <- strsplit(
        ifelse(is.na(variant_analysis$significant_cell_types) | variant_analysis$significant_cell_types == "",
          "", variant_analysis$significant_cell_types
        ),
        ",",
        fixed = TRUE
      )
      feat_cts_split <- strsplit(
        ifelse(is.na(variant_analysis$significant_cts_str) | variant_analysis$significant_cts_str == "",
          "", variant_analysis$significant_cts_str
        ),
        ",",
        fixed = TRUE
      )

      # Determine if each variant is shared (vectorized)
      variant_analysis[, `:=`(
        # Count overlapping cell types
        overlap_count = mapply(function(var_cts, feat_cts) {
          if (length(var_cts) == 0L || (length(var_cts) == 1L && var_cts[1L] == "")) {
            return(0L)
          }
          if (length(feat_cts) == 0L || (length(feat_cts) == 1L && feat_cts[1L] == "")) {
            return(0L)
          }
          length(intersect(var_cts, feat_cts))
        }, var_cts_split, feat_cts_split),

        # Check if feature or variant is underpowered
        feature_is_underpowered = specificity_category == "Likely shared but underpowered",
        variant_is_underpowered = variant_specificity == "Likely shared but underpowered",

        # Check if feature is truly shared (not underpowered)
        feature_is_shared = specificity_category == "Cross-lineage shared"
      )]

      # Determine if shared based on hierarchy
      # Logic:
      # 1. If overlap < 2, not shared
      # 2. If feature is underpowered, be lenient (any overlapping variant is shared)
      # 3. If variant is underpowered AND feature is truly shared, also be lenient
      # 4. Otherwise, require exact specificity match
      variant_analysis[, is_shared := ifelse(
        overlap_count >= 2,
        ifelse(
          feature_is_underpowered,
          TRUE, # Feature underpowered: be lenient
          ifelse(
            variant_is_underpowered & feature_is_shared,
            TRUE, # Variant underpowered but feature is truly shared: be lenient
            variant_specificity == specificity_category # Otherwise: require exact match
          )
        ),
        FALSE
      )]

      # Aggregate patterns by feature
      feature_patterns <- variant_analysis[,
        {
          if (.N == 0 || all(is.na(cochran_q_pval))) {
            # No variants or no heterogeneity tests available - return NA
            list(pattern = NA_character_)
          } else if (!any(is_shared)) {
            # All variants are distinct (not shared with feature)
            list(pattern = "distinct_variants")
          } else {
            # Get shared variants
            shared_data <- .SD[is_shared == TRUE]
            if (nrow(shared_data) == 0) {
              list(pattern = "distinct_variants")
            } else {
              has_opposite <- any((shared_data$cochran_q_pval < cochran_q_threshold) & (shared_data$direction == "+-"), na.rm = TRUE)
              has_heterogeneous <- any(shared_data$cochran_q_pval < cochran_q_threshold, na.rm = TRUE)

              if (has_opposite) {
                list(pattern = "shared_opposite")
              } else if (has_heterogeneous) {
                list(pattern = "shared_heterogeneous")
              } else {
                list(pattern = "shared_consistent")
              }
            }
          }
        },
        by = feature_id
      ]

      # Initialize all patterns as NA
      results$hierarchical_variant_pattern <- NA_character_

      # Set patterns for features with no variants
      no_variant_features <- features_to_process[!feature_id %in% meta_filtered$feature_id]
      if (nrow(no_variant_features) > 0) {
        results$hierarchical_variant_pattern[no_variant_features$idx] <- NA_character_
      }

      # Set patterns from analysis (vectorized join instead of row-by-row loop)
      pattern_map <- merge(
        feature_patterns[, .(feature_id, pattern)],
        features_dt[, .(feature_id, idx)],
        by = "feature_id"
      )
      if (nrow(pattern_map) > 0) {
        results$hierarchical_variant_pattern[pattern_map$idx] <- pattern_map$pattern
      }
    } else {
      results$hierarchical_variant_pattern <- NA_character_
    }

    # Print summary of hierarchical patterns
    hier_pattern_counts <- table(results$hierarchical_variant_pattern, useNA = "ifany")
    cli::cli_alert_info("Hierarchical pattern distribution:")
    for (p in names(hier_pattern_counts)) {
      if (!is.na(p)) {
        cli::cli_bullets(setNames(paste0(p, ": ", hier_pattern_counts[[p]]), " "))
      }
    }
  }

  # Extract variant details from SuSiE results
  variant_details <- NULL

  # We need to extract SuSiE data for variants used in heterogeneity analysis
  # This requires the SuSiE results which should be passed to this function
  if (!is.null(susie_results) && nrow(meta_filtered) > 0) {
    cli::cli_alert_info("Extracting SuSiE variant details for heterogeneity analysis")

    # Create variant-level patterns from meta_filtered
    # Each variant has its own pattern based on its Cochran's Q p-value and direction
    variant_patterns <- meta_filtered[, .(
      variant_id,
      feature_id,
      cochran_q_pval, # Keep the p-value for output
      variant_pattern = ifelse(
        is.na(cochran_q_pval), "shared_consistent",
        ifelse(
          cochran_q_pval < cochran_q_threshold,
          ifelse(direction == "+-", "shared_opposite", "shared_heterogeneous"),
          "shared_consistent"
        )
      )
    )]

    # Get feature-level patterns and top variants
    feature_patterns <- feature_dt[!is.na(pattern), .(
      feature_id,
      feature_pattern = pattern,
      top_variant,
      max_pip,
      max_chisq
    )]

    # Also get hierarchical patterns if they exist
    hierarchical_patterns_dt <- NULL
    if ("hierarchical_variant_pattern" %in% names(results) && !all(is.na(results$hierarchical_variant_pattern))) {
      hierarchical_patterns_dt <- data.table::data.table(
        feature_id = feature_ids,
        hierarchical_pattern = results$hierarchical_variant_pattern
      )
    }

    # Create a combined lookup table with both feature and variant patterns
    pattern_lookup <- merge(
      variant_patterns,
      feature_patterns,
      by = "feature_id",
      all.x = TRUE
    )

    # Add hierarchical patterns if available
    if (!is.null(hierarchical_patterns_dt)) {
      pattern_lookup <- merge(
        pattern_lookup,
        hierarchical_patterns_dt,
        by = "feature_id",
        all.x = TRUE
      )
    } else {
      pattern_lookup[, hierarchical_pattern := NA_character_]
    }

    pattern_lookup[, is_top_variant := variant_id == top_variant]

    # Extract SuSiE data using vectorized operations
    variant_details_list <- list()

    # Process all cell types at once
    for (ct in names(susie_results)) {
      ct_results <- susie_results[[ct]]

      # Get all features that have data in this cell type
      available_features <- names(ct_results)
      features_to_extract <- unique(meta_filtered$feature_id[meta_filtered$feature_id %in% available_features])

      if (length(features_to_extract) == 0) next

      # Extract all SuSiE data for this cell type in one go
      ct_data_list <- lapply(features_to_extract, function(feat_id) {
        feat_data <- ct_results[[feat_id]]
        if (is.null(feat_data$variant_names)) {
          return(NULL)
        }

        # Find variants to extract (those with max_pip > 0.5)
        feat_variants <- meta_filtered[feature_id == feat_id, unique(variant_id)]
        variant_indices <- which(feat_data$variant_names %in% feat_variants)

        if (length(variant_indices) == 0) {
          return(NULL)
        }

        # Create data.table with SuSiE data
        data.table::data.table(
          gene_id = feat_id,
          cell_type = ct,
          variant_id = feat_data$variant_names[variant_indices],
          pip = if (!is.null(feat_data$pip)) feat_data$pip[variant_indices] else NA_real_,
          beta = if (!is.null(feat_data$beta)) feat_data$beta[variant_indices] else NA_real_,
          se = if (!is.null(feat_data$se)) feat_data$se[variant_indices] else NA_real_
        )
      })

      # Combine all features for this cell type
      ct_data <- data.table::rbindlist(ct_data_list[!sapply(ct_data_list, is.null)])

      if (nrow(ct_data) > 0) {
        variant_details_list[[length(variant_details_list) + 1]] <- ct_data
      }
    }

    # Combine all cell types
    if (length(variant_details_list) > 0) {
      variant_details <- data.table::rbindlist(variant_details_list, fill = TRUE)

      # Join with pattern lookup to add both feature and variant patterns
      # Use gene_id as the join key for feature_id
      pattern_lookup_for_join <- copy(pattern_lookup)
      data.table::setnames(pattern_lookup_for_join, "feature_id", "gene_id")

      # Set keys for efficient joining
      data.table::setkey(variant_details, gene_id, variant_id)
      data.table::setkey(pattern_lookup_for_join, gene_id, variant_id)

      # Merge patterns and top variant information
      variant_details <- merge(
        variant_details,
        pattern_lookup_for_join[, .(gene_id, variant_id, feature_pattern, variant_pattern, hierarchical_pattern, cochran_q_pval, is_top_variant)],
        by = c("gene_id", "variant_id"),
        all.x = TRUE
      )

      # Add gene specificity from results
      gene_specificity_dt <- data.table::data.table(
        feature_id = feature_ids,
        gene_specificity = spec_values
      )
      # Rename column for joining
      if (data_type == "gene") {
        data.table::setnames(gene_specificity_dt, "feature_id", "gene_id")
      } else {
        data.table::setnames(gene_specificity_dt, "feature_id", "peak_id")
        data.table::setnames(variant_details, "gene_id", "peak_id")
      }

      # Merge gene/feature specificity
      variant_details <- merge(
        variant_details,
        gene_specificity_dt,
        by = if (data_type == "gene") "gene_id" else "peak_id",
        all.x = TRUE
      )

      # Add variant specificity from variant_feature_specificity
      if (!is.null(variant_feature_specificity) && nrow(variant_feature_specificity) > 0) {
        # Get unique variant specificity for each variant-feature pair
        vfs_for_merge <- variant_feature_specificity[
          feature_type == data_type,
          .(variant_specificity = first(variant_specificity)),
          by = .(feature_id, variant_id)
        ]

        # Rename feature_id column to match
        if (data_type == "gene") {
          data.table::setnames(vfs_for_merge, "feature_id", "gene_id")
        } else {
          data.table::setnames(vfs_for_merge, "feature_id", "peak_id")
        }

        # Merge variant specificity
        variant_details <- merge(
          variant_details,
          vfs_for_merge,
          by = c(if (data_type == "gene") "gene_id" else "peak_id", "variant_id"),
          all.x = TRUE
        )
      } else {
        variant_details[, variant_specificity := NA_character_]
      }

      # Handle missing patterns (shouldn't happen but just in case)
      variant_details[is.na(feature_pattern), feature_pattern := "unknown"]
      variant_details[is.na(variant_pattern), variant_pattern := "unknown"]
      variant_details[is.na(hierarchical_pattern), hierarchical_pattern := "unknown"]
      variant_details[is.na(is_top_variant), is_top_variant := FALSE]
      variant_details[is.na(gene_specificity), gene_specificity := "unknown"]
      variant_details[is.na(variant_specificity), variant_specificity := "unknown"]

      # Rename back if peak
      if (data_type == "peak") {
        data.table::setnames(variant_details, "peak_id", "gene_id")
        data.table::setnames(variant_details, "gene_specificity", "peak_specificity")
      }

      # Reorder columns as requested
      if (data_type == "gene") {
        variant_details <- variant_details[, .(
          gene_id,
          cell_type,
          variant_id,
          pip,
          beta,
          se,
          cochran_q_pval, # Cochran's Q p-value for heterogeneity
          gene_specificity, # Gene specificity category
          variant_specificity, # Variant specificity category
          feature_pattern, # Feature-level heterogeneity pattern (variant_heterogeneity_pattern)
          hierarchical_pattern, # Hierarchical variant pattern
          variant_pattern, # Variant-level heterogeneity pattern
          is_top_variant
        )]
      } else {
        variant_details <- variant_details[, .(
          gene_id,
          cell_type,
          variant_id,
          pip,
          beta,
          se,
          cochran_q_pval, # Cochran's Q p-value for heterogeneity
          peak_specificity, # Peak specificity category
          variant_specificity, # Variant specificity category
          feature_pattern, # Feature-level heterogeneity pattern (variant_heterogeneity_pattern)
          hierarchical_pattern, # Hierarchical variant pattern
          variant_pattern, # Variant-level heterogeneity pattern
          is_top_variant
        )]
      }

      # Sort by gene_id, variant_id, and cell_type for consistent output
      data.table::setorder(variant_details, gene_id, variant_id, cell_type)
    }
  }

  return(list(
    results = results,
    variant_details = variant_details
  ))
}

#' Extract CS Cluster Variant Details
#'
#' Creates a detailed data table of top variants for each CS cluster
#'
#' @param cluster_top_variants Data table with top variant per cluster from analyze_cs_cluster_heterogeneity
#' @param cs_clusters_filtered Data table with filtered CS cluster mappings
#' @param feature_type Either "gene" or "peak"
#' @return Data table with detailed variant information per cluster
#' @keywords internal
extract_cs_cluster_variant_details <- function(cluster_top_variants, cs_clusters_filtered,
                                               feature_type = "gene") {
  # If no cluster top variants, return empty data.table
  if (is.null(cluster_top_variants) || nrow(cluster_top_variants) == 0) {
    return(data.table::data.table())
  }

  # Get cell types per cluster
  cluster_cell_types <- cs_clusters_filtered[, .(
    cell_types = paste(unique(sort(cell_type)), collapse = ","),
    n_cell_types = length(unique(cell_type))
  ), by = .(feature_id, cluster_id)]

  # Join cluster top variants with cell type information
  variant_details <- merge(
    cluster_top_variants,
    cluster_cell_types,
    by.x = c("feature_ids", "cluster_id"),
    by.y = c("feature_id", "cluster_id"),
    all.x = TRUE
  )

  # Rename columns appropriately
  id_col <- if (feature_type == "gene") "gene_id" else "peak_id"
  data.table::setnames(variant_details, "feature_ids", id_col)

  # Add feature-level pattern (to be joined from main results later)
  variant_details[, cluster_pattern := pattern]

  # Select and order columns for output
  output_cols <- c(
    id_col,
    "cluster_id",
    "variant_id",
    "max_chisq",
    "max_pip",
    "cochran_q_pval",
    "direction",
    "cluster_pattern",
    "cell_types",
    "n_cell_types"
  )

  # Keep only columns that exist
  output_cols <- intersect(output_cols, names(variant_details))
  variant_details <- variant_details[, ..output_cols]

  # Sort by feature, cluster
  data.table::setorderv(variant_details, c(id_col, "cluster_id"))

  return(variant_details)
}

#' Extract Per-Cell-Type SuSiE Results for CS Cluster Variants
#'
#' Extracts detailed per-cell-type SuSiE results (beta, se, pip) for all variants
#' in CS clusters using vectorized operations
#'
#' @param features Vector of feature IDs to extract
#' @param cs_clusters Data table with CS to cluster mappings
#' @param cs_cluster_variants Data table with cluster to variant mappings
#' @param susie_results List of SuSiE results by cell type
#' @param feature_type Either "gene" or "peak"
#' @return Data table with per-cell-type SuSiE results for CS cluster variants
extract_cs_cluster_susie_details <- function(features, cs_clusters, cs_cluster_variants,
                                             susie_results, feature_type = "gene") {
  # Input validation
  if (is.null(features) || length(features) == 0 ||
    is.null(susie_results) || length(susie_results) == 0) {
    return(data.table::data.table())
  }

  # Ensure data.tables
  if (!data.table::is.data.table(cs_clusters)) data.table::setDT(cs_clusters)
  if (!data.table::is.data.table(cs_cluster_variants)) data.table::setDT(cs_cluster_variants)

  # Filter to features of interest
  cs_clusters_filtered <- cs_clusters[feature_id %in% features]
  if (nrow(cs_clusters_filtered) == 0) {
    return(data.table::data.table())
  }

  # Get all relevant clusters
  all_clusters <- unique(cs_clusters_filtered$cluster_id)

  # Get variants for these clusters and features
  cluster_vars_filtered <- cs_cluster_variants[cluster_id %in% all_clusters & feature_ids %in% features]
  if (nrow(cluster_vars_filtered) == 0) {
    return(data.table::data.table())
  }

  # Pre-compute target variant set and feature set for fast filtering
  target_variants <- unique(cluster_vars_filtered$variant_id)
  target_features <- unique(cluster_vars_filtered$feature_ids)

  # Prepare cluster lookup table (keyed for fast joins)
  cluster_lookup <- unique(cluster_vars_filtered[, .(feature_ids, variant_id, cluster_id)])
  data.table::setkey(cluster_lookup, feature_ids, variant_id)

  # Flatten SuSiE data into one data.table per cell type.
  # Per-feature %chin% filtering has per-call overhead with millions of target
  # variants, so instead we:
  # 1. Extract ALL variant_names from ALL features in one pass (lapply + unlist)
  # 2. Do a SINGLE %chin% filter on the concatenated vector
  # 3. Use index-based subsetting to reconstruct the filtered data
  # This reduces R-level overhead from O(n_features) to O(1) per cell type.

  # Process each cell type (parallel if mc.cores > 1 on Unix)
  n_cores <- getOption("mc.cores", 1L)
  use_parallel <- .Platform$OS.type == "unix" && n_cores > 1 && length(susie_results) > 1

  process_ct <- function(ct) {
    ct_results <- susie_results[[ct]]

    # Get features that have results in this cell type AND are in our target set
    available_features <- intersect(names(ct_results), target_features)

    if (length(available_features) == 0) {
      return(NULL)
    }

    # Step 1: Extract ALL variant_names across features (one lapply, no per-feature overhead)
    all_vnames <- lapply(ct_results[available_features], `[[`, "variant_names")

    # Compute lengths for rep() and subsetting
    lens <- lengths(all_vnames)
    non_empty <- lens > 0
    if (!any(non_empty)) {
      return(NULL)
    }

    # Filter to non-empty only
    available_features <- available_features[non_empty]
    all_vnames <- all_vnames[non_empty]
    lens <- lens[non_empty]

    # Step 2: Concatenate all variant names into one vector
    all_vn <- unlist(all_vnames, use.names = FALSE)

    # Step 3: Single %chin% filter on concatenated vector
    keep <- all_vn %chin% target_variants
    if (!any(keep)) {
      return(NULL)
    }

    # Step 4: Build feature_id vector via rep
    feature_ids <- rep(available_features, lens)

    # Step 5: Extract beta/se/pip in bulk (same structure)
    all_beta <- unlist(lapply(ct_results[available_features], function(fd) {
      if (!is.null(fd$beta)) fd$beta else rep(NA_real_, length(fd$variant_names))
    }), use.names = FALSE)
    all_se <- unlist(lapply(ct_results[available_features], function(fd) {
      if (!is.null(fd$se)) fd$se else rep(NA_real_, length(fd$variant_names))
    }), use.names = FALSE)
    all_pip <- unlist(lapply(ct_results[available_features], function(fd) {
      if (!is.null(fd$pip)) fd$pip else rep(NA_real_, length(fd$variant_names))
    }), use.names = FALSE)

    # Step 6: Apply single filter to all vectors at once
    ct_data <- data.table::data.table(
      variant_id = all_vn[keep],
      beta = all_beta[keep],
      se = all_se[keep],
      pip = all_pip[keep],
      feature_id = feature_ids[keep]
    )

    # Single keyed join to add cluster info
    ct_data <- cluster_lookup[ct_data,
      on = .(feature_ids = feature_id, variant_id),
      nomatch = NULL,
      .(
        feature_id = feature_ids, variant_id, cluster_id, beta = i.beta,
        se = i.se, pip = i.pip
      )
    ]

    if (nrow(ct_data) == 0) {
      return(NULL)
    }

    # Add cell type and compute derived columns vectorized
    ct_data[, cell_type := ct]
    ct_data[, `:=`(
      z_score = fifelse(se > 0, beta / se, NA_real_),
      p_value = fifelse(se > 0, 2 * pnorm(-abs(beta / se)), NA_real_)
    )]

    return(ct_data)
  }

  if (use_parallel) {
    susie_data_list <- parallel::mclapply(names(susie_results), process_ct,
      mc.cores = min(n_cores, length(susie_results))
    )
  } else {
    susie_data_list <- lapply(names(susie_results), process_ct)
  }

  # Combine all cell types
  susie_details <- data.table::rbindlist(
    susie_data_list[!vapply(susie_data_list, is.null, logical(1))],
    use.names = TRUE
  )

  if (nrow(susie_details) == 0) {
    return(data.table::data.table())
  }

  # Rename feature_id column appropriately
  id_col <- if (feature_type == "gene") "gene_id" else "peak_id"
  data.table::setnames(susie_details, "feature_id", id_col)

  # Reorder columns for output
  col_order <- c(
    id_col, "cluster_id", "variant_id", "cell_type",
    "beta", "se", "pip", "z_score", "p_value"
  )
  data.table::setcolorder(susie_details, intersect(col_order, names(susie_details)))

  # Sort by feature, cluster, variant, cell type
  data.table::setorderv(susie_details, c(id_col, "cluster_id", "variant_id", "cell_type"))

  return(susie_details)
}

#' Analyze CS Cluster-based Heterogeneity
#'
#' Analyzes heterogeneity patterns using CS cluster information for more accurate
#' detection of shared vs distinct signals across cell types
#'
#' @param features Vector of feature IDs to analyze
#' @param cs_clusters Data table with CS to cluster mappings
#' @param cs_cluster_variants Data table with cluster to variant mappings
#' @param meta_data Data table with variant heterogeneity information
#' @param cochran_q_threshold P-value threshold for Cochran's Q test
#' @param feature_type Either "gene" or "peak"
#' @param return_variant_details If TRUE, returns list with results and variant details; otherwise just results
#' @return Data table with heterogeneity analysis results or list with results and variant details
#' @keywords internal
analyze_cs_cluster_heterogeneity <- function(features, cs_clusters, cs_cluster_variants,
                                             meta_data, cochran_q_threshold = 5e-8,
                                             feature_type = "gene",
                                             return_variant_details = FALSE) {
  # Input validation
  if (is.null(features) || length(features) == 0) {
    return(data.table::data.table())
  }

  if (is.null(cs_clusters) || nrow(cs_clusters) == 0 ||
    is.null(cs_cluster_variants) || nrow(cs_cluster_variants) == 0 ||
    is.null(meta_data) || nrow(meta_data) == 0) {
    # Return empty results if any required data is missing
    n_features <- length(features)
    return(data.table::data.table(
      feature_id = features,
      variant_heterogeneity_pattern = rep(NA_character_, n_features),
      n_cs_clusters = rep(0L, n_features),
      top_variant = rep(NA_character_, n_features),
      max_chisq = rep(NA_real_, n_features),
      variant_heterogeneity = rep(FALSE, n_features)
    ))
  }

  # Ensure all inputs are data.tables with proper keys
  if (!data.table::is.data.table(cs_clusters)) data.table::setDT(cs_clusters)
  if (!data.table::is.data.table(cs_cluster_variants)) data.table::setDT(cs_cluster_variants)
  if (!data.table::is.data.table(meta_data)) data.table::setDT(meta_data)

  # Filter to only features we're analyzing
  cs_clusters_filtered <- cs_clusters[feature_id %in% features]

  if (nrow(cs_clusters_filtered) == 0) {
    # No clusters for any of the features
    return(data.table::data.table(
      feature_id = features,
      variant_heterogeneity_pattern = rep(NA_character_, length(features)),
      n_cs_clusters = rep(0L, length(features)),
      top_variant = rep(NA_character_, length(features)),
      max_chisq = rep(NA_real_, length(features)),
      variant_heterogeneity = rep(FALSE, length(features))
    ))
  }

  # Step 1: Count clusters per feature and get unique clusters
  feature_cluster_counts <- cs_clusters_filtered[, .(
    n_cs_clusters = length(unique(cluster_id)),
    clusters = list(unique(cluster_id))
  ), by = feature_id]

  # Step 2: Join cluster variants with meta data to get top variant per cluster
  # First, get all relevant clusters
  all_clusters <- unique(cs_clusters_filtered$cluster_id)

  # Get variants for these clusters AND filter to only features we're analyzing
  # This is crucial to avoid cross-contamination between features
  cluster_vars_filtered <- cs_cluster_variants[cluster_id %in% all_clusters & feature_ids %in% features]

  # Join with meta data using keyed binary merge for speed
  # Set keys on both sides so data.table uses fast binary join
  if (!identical(data.table::key(meta_data), c("variant_id", "feature_id"))) {
    data.table::setkey(meta_data, variant_id, feature_id)
  }
  data.table::setkey(cluster_vars_filtered, variant_id, feature_ids)
  cluster_meta <- meta_data[cluster_vars_filtered,
    on = .(variant_id = variant_id, feature_id = feature_ids),
    nomatch = NULL
  ]
  # Restore feature_ids column name for downstream cs_cluster_variants compatibility
  data.table::setnames(cluster_meta, "feature_id", "feature_ids")

  if (nrow(cluster_meta) == 0) {
    # No matching meta data for any clusters
    results_dt <- data.table::data.table(feature_id = features)
    results_dt <- merge(
      results_dt,
      feature_cluster_counts,
      by = "feature_id",
      all.x = TRUE
    )
    results_dt[is.na(n_cs_clusters), n_cs_clusters := 0L]
    results_dt[, `:=`(
      variant_heterogeneity_pattern = NA_character_,
      top_variant = NA_character_,
      max_chisq = NA_real_,
      variant_heterogeneity = FALSE,
      clusters = NULL
    )]
    return(results_dt)
  }

  # Step 3: Find top variant per cluster based on max_chisq
  # Use .I[which.max(...)] instead of .SD[which.max(...)] to avoid materializing
  # a sub-data.table per group (294K groups -> 38s with .SD vs <1s with .I)
  cluster_meta_valid <- cluster_meta[!is.na(max_chisq)]
  top_idx <- cluster_meta_valid[, .I[which.max(max_chisq)], by = .(cluster_id, feature_ids)]$V1
  cluster_top_variants <- cluster_meta_valid[top_idx]

  # Get number of cell types per cluster to determine if heterogeneity is applicable
  cluster_n_cts <- cs_clusters_filtered[, .(n_cell_types = length(unique(cell_type))), by = .(feature_id, cluster_id)]
  cluster_top_variants <- merge(cluster_top_variants, cluster_n_cts,
    by.x = c("feature_ids", "cluster_id"),
    by.y = c("feature_id", "cluster_id"),
    all.x = TRUE
  )

  # Determine pattern for each cluster
  # For single cell type clusters, heterogeneity is not applicable (NA)
  cluster_top_variants[, pattern := data.table::fcase(
    n_cell_types == 1, NA_character_, # Single cell type - heterogeneity not applicable
    !is.na(cochran_q_pval) & cochran_q_pval < cochran_q_threshold & !is.na(direction) & direction %in% c("+-", "-+"), "opposite",
    !is.na(cochran_q_pval) & cochran_q_pval < cochran_q_threshold, "heterogeneous",
    default = "consistent"
  )]

  # Step 4: Aggregate to feature level
  # Split into two fast operations instead of one .SD[which.max()] per group

  # 4a: Get top variant per feature using .I indexing (avoids .SD materialization)
  top_feat_idx <- cluster_top_variants[, .I[which.max(max_chisq)], by = feature_ids]$V1
  feature_top <- cluster_top_variants[top_feat_idx, .(feature_ids, top_variant = variant_id, max_chisq)]

  # 4b: Compute pattern per feature using vectorized priority logic
  feature_pattern_dt <- cluster_top_variants[!is.na(pattern), .(
    has_opposite = any(pattern == "opposite"),
    has_heterogeneous = any(pattern == "heterogeneous")
  ), by = feature_ids]
  feature_pattern_dt[, variant_heterogeneity_pattern := data.table::fcase(
    has_opposite, "shared_opposite",
    has_heterogeneous, "shared_heterogeneous",
    default = "shared_consistent"
  )]

  # 4c: Merge top variants with patterns
  feature_patterns <- merge(feature_top, feature_pattern_dt[, .(feature_ids, variant_heterogeneity_pattern)],
    by = "feature_ids", all.x = TRUE
  )
  # Features with only single-CT clusters have no pattern row — leave as NA
  feature_patterns[is.na(variant_heterogeneity_pattern), variant_heterogeneity_pattern := NA_character_]

  # Step 5: Check for distinct variants (multiple clusters with no cell type overlap)
  # Get cell type counts per cluster for multi-cluster features
  multi_cluster_features <- feature_cluster_counts[n_cs_clusters > 1, feature_id]

  if (length(multi_cluster_features) > 0) {
    cluster_cell_types <- cs_clusters_filtered[feature_id %in% multi_cluster_features, .(
      n_cell_types = length(unique(cell_type))
    ), by = .(feature_id, cluster_id)]

    # Check if any feature has clusters that don't share cell types
    distinct_features <- cluster_cell_types[, .(
      max_cell_types = max(n_cell_types)
    ), by = feature_id][max_cell_types == 1, feature_id]

    # Update pattern for these features
    if (length(distinct_features) > 0) {
      feature_patterns[
        feature_ids %in% distinct_features,
        variant_heterogeneity_pattern := "distinct_variants"
      ]
    }
  }

  # Step 6: Combine all results
  results_dt <- data.table::data.table(feature_id = features)

  # Merge with cluster counts
  results_dt <- merge(results_dt, feature_cluster_counts[, .(feature_id, n_cs_clusters)],
    by = "feature_id", all.x = TRUE
  )

  # Merge with patterns
  results_dt <- merge(results_dt, feature_patterns,
    by.x = "feature_id", by.y = "feature_ids", all.x = TRUE
  )

  # Fill in missing values
  results_dt[is.na(n_cs_clusters), n_cs_clusters := 0L]
  results_dt[is.na(variant_heterogeneity_pattern), variant_heterogeneity_pattern :=
    ifelse(n_cs_clusters > 0, "distinct_variants", NA_character_)]
  results_dt[is.na(top_variant), top_variant := NA_character_]
  results_dt[is.na(max_chisq), max_chisq := NA_real_]

  # Add heterogeneity flag
  results_dt[, variant_heterogeneity := variant_heterogeneity_pattern %in%
    c("shared_opposite", "shared_heterogeneous")]

  # If variant details are requested, extract them
  if (return_variant_details) {
    # Extract detailed variant information for each cluster
    variant_details <- extract_cs_cluster_variant_details(
      cluster_top_variants = cluster_top_variants,
      cs_clusters_filtered = cs_clusters_filtered,
      feature_type = feature_type
    )

    return(list(
      results = results_dt,
      variant_details = variant_details
    ))
  }

  return(results_dt)
}

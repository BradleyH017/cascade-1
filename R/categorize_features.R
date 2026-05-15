#' Categorize features (genes or peaks)
#'
#' High-level dispatch for categorizing genes or peaks across cell types,
#' powered by vectorized C++ kernels.
#'
#' @param feature_data Feature data from load_feature_data(config, feature_type, chromosomes)
#' @param lfsr_results Pre-computed LFSR results from load_lfsr_results()
#' @param susie_results SuSiE results (optional)
#' @param meta_data Required. Pre-computed meta data with Cochran's Q values
#' @param variant_feature_specificity Required. Variant-feature specificity mapping from variant categorization
#' @param feature_type Either "gene" or "peak"
#' @param hierarchy A CellTypeHierarchy object; defaults to DEFAULT_CELL_HIERARCHY
#' @param lfsr_sig_threshold Significance threshold used for both ACAT q-value
#'   significance and LFSR gray zone lower bound (by design, both use the same
#'   threshold). This is NOT acat_fdr_threshold (which is used only in variant
#'   data loading).
#' @param lfsr_null_threshold LFSR null hypothesis threshold
#' @param cochran_q_threshold P-value threshold for Cochran's Q heterogeneity test (default 5e-8)
#' @param use_cs_clusters If TRUE, use credible-set cluster assignments when available; otherwise rely on Cochran's Q only
#' @param cs_clusters Optional credible-set cluster assignments (loaded by load_cs_clusters())
#' @param cs_cluster_variants Optional cluster-level variant details (loaded by load_cs_cluster_variants())
#' @return List with two elements: categories (data frame with categorization results) and variant_details (data frame with variant details)
#' @export
categorize_features <- function(feature_data,
                                lfsr_results,
                                susie_results,
                                meta_data,
                                variant_feature_specificity,
                                feature_type = "gene",
                                hierarchy = DEFAULT_CELL_HIERARCHY,
                                lfsr_sig_threshold = LFSR_SIG_THRESHOLD,
                                lfsr_null_threshold = LFSR_NULL_THRESHOLD,
                                cochran_q_threshold = 5e-8,
                                use_cs_clusters = TRUE,
                                cs_clusters = NULL,
                                cs_cluster_variants = NULL) {
  cli::cli_alert_info("Starting {feature_type} categorization")


  # Use the unified approach
  cli::cli_alert_info("Using unified categorization approach")

  # Use pre-computed feature-level LFSR lookup from load_lfsr_results
  lfsr_lookup <- NULL
  if (!is.null(lfsr_results)) {
    lfsr_key <- if (feature_type == "gene") "eqtl_feature" else "caqtl_feature"
    lfsr_lookup <- lfsr_results[[lfsr_key]]
  }

  cpp_params <- hierarchy_to_cpp_params(hierarchy)
  results <- categorize_features_from_acat(
    acat_matrix = feature_data$acat_matrix,
    lfsr_lookup = lfsr_lookup,
    feature_type = feature_type,
    l2_to_l1_mapping = cpp_params$l2_to_l1_mapping,
    lineage_groups = cpp_params$lineage_groups,
    subgroup_levels = cpp_params$subgroup_levels,
    bulk_cts = cpp_params$bulk_cts,
    other_cts = cpp_params$other_cts,
    sig_threshold = lfsr_sig_threshold,
    lfsr_null_threshold = lfsr_null_threshold,
    specificity_categories = cpp_params$specificity_categories
  )

  # Rename columns to match expected output
  names(results)[names(results) == "primary_category"] <- "cell_type_specificity"

  # Rename feature_id to appropriate column name
  id_col <- paste0(feature_type, "_id")
  names(results)[names(results) == "feature_id"] <- id_col

  # Add variant analysis (SuSiE results are now guaranteed from load_feature_data)
  cli::cli_alert_info("Adding variant pattern analysis")
  # Always use Cochran's Q method with extended extraction (or CS clusters if available)
  variant_result <- add_variant_analysis(results, susie_results, feature_type, lfsr_results,
    meta_data, variant_feature_specificity,
    cochran_q_threshold = cochran_q_threshold,
    use_cs_clusters = use_cs_clusters,
    cs_clusters = cs_clusters,
    cs_cluster_variants = cs_cluster_variants
  )

  # No need to handle NULL anymore since add_variant_analysis always returns a list

  cli::cli_alert_success("Categorized {nrow(variant_result$results)} {feature_type}s")

  # Return both results and variant details as a list
  return(list(
    categories = variant_result$results,
    variant_details = variant_result$variant_details
  ))
}

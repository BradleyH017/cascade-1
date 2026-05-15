#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom stats pnorm setNames
#' @importFrom utils head modifyList packageVersion write.table
#' @useDynLib cascade, .registration = TRUE
## usethis namespace: end
NULL


# Silence R CMD check NOTEs about data.table NSE / non-standard column names
# referenced inside data.table[] expressions.
utils::globalVariables(c(
  ".",
  "..output_cols", "..select_cols",
  "affected_cell_types", "associated_genes", "associated_peaks",
  "best_cell_types", "caqtl", "caqtl_peaks", "cascade_peak_genes",
  "cell_type", "cell_type_specificity", "cluster_id", "cluster_pattern",
  "cochran_q_pval", "config_cell_type", "ct", "ct_idx", "direction",
  "eqtl", "eqtl_genes", "feature", "feature_id", "feature_ids",
  "feature_is_shared", "feature_is_underpowered", "feature_pattern",
  "feature_type", "features", "features_list",
  "gene_affected_cell_types", "gene_cell_type_specificity", "gene_id",
  "gene_specificity", "has_heterogeneous", "has_opposite",
  "hierarchical_pattern",
  "i.beta", "i.lfsr", "i.lfsr_str", "i.lfsr_value", "i.pip", "i.se",
  "idx", "is_shared", "is_top_variant", "lfsr", "lfsr_value",
  "link_mechanism", "list_idx", "max_cell_types", "max_chisq", "max_pip",
  "mechanism", "min_lfsr", "min_lfsr_gene", "min_lfsr_peak",
  "n_affected_cell_types", "n_caqtl_peaks", "n_cell_types",
  "n_cs_clusters", "n_eqtl_genes", "n_gene_affected_cell_types",
  "n_peak_affected_cell_types", "n_sig_cts", "n_significant_cts",
  "n_variants", "orig_idx", "orig_row", "overlap_count", "pair", "pairs",
  "pairs_list", "pattern", "peak_affected_cell_types",
  "peak_cell_type_specificity", "peak_gene_link", "peak_overlap",
  "peak_specificity", "pip", "process",
  "qtl_mechanism_category", "qtl_pattern", "qtl_pattern_number",
  "qtl_type", "se", "significant_cell_types", "significant_cts",
  "significant_cts_str", "specificity_category", "top_variant",
  "var_idx", "variant_heterogeneity", "variant_heterogeneity_pattern",
  "variant_id", "variant_is_underpowered", "variant_pattern",
  "variant_specificity"
))

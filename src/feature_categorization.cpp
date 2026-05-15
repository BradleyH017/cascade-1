#include <Rcpp.h>

#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "cascade_utils.h"

using namespace Rcpp;

// [[Rcpp::export]]
std::string categorize_cell_specificity(
    CharacterVector significant_cts, CharacterVector specificity_categories,
    Nullable<NumericVector> lfsr_values = R_NilValue,
    Nullable<CharacterVector> lfsr_names = R_NilValue,
    List lineage_groups = List::create(), List subgroup_levels = List::create(),
    CharacterVector bulk_cts = CharacterVector::create(),
    CharacterVector other_cts = CharacterVector::create(),
    double lfsr_sig_threshold = 0.05, double lfsr_null_threshold = 0.5) {
  // Build generalized cell type sets
  cascade::CellTypeSets ct_sets(lineage_groups, subgroup_levels, bulk_cts,
                                other_cts);

  // Build LFSR map if provided
  std::unordered_map<std::string, double> lfsr_map;
  bool has_lfsr = lfsr_values.isNotNull() && lfsr_names.isNotNull();
  if (has_lfsr) {
    NumericVector lfsr_vals = as<NumericVector>(lfsr_values);
    CharacterVector lfsr_nms = as<CharacterVector>(lfsr_names);

    if (lfsr_vals.size() != lfsr_nms.size()) {
      stop("LFSR values and names must have the same length");
    }

    for (int i = 0; i < lfsr_vals.size(); i++) {
      if (!NumericVector::is_na(lfsr_vals[i])) {
        lfsr_map[as<std::string>(lfsr_nms[i])] = lfsr_vals[i];
      }
    }
  }

  // Convert significant cell types to vector
  std::vector<std::string> sig_cts = cascade::to_string_vector(significant_cts);

  // Convert specificity categories to vector
  std::vector<std::string> categories =
      cascade::to_string_vector(specificity_categories);

  // Call generalized internal function
  return cascade::categorize_single_feature_internal(
      sig_cts, lfsr_map, ct_sets, lfsr_sig_threshold, lfsr_null_threshold,
      categories);
}

// [[Rcpp::export]]
DataFrame categorize_features_from_acat(
    DataFrame acat_matrix, CharacterVector specificity_categories,
    Nullable<DataFrame> lfsr_lookup = R_NilValue,
    std::string feature_type = "gene",
    Nullable<List> l2_to_l1_mapping = R_NilValue,
    List lineage_groups = List::create(), List subgroup_levels = List::create(),
    CharacterVector bulk_cts = CharacterVector::create(),
    CharacterVector other_cts = CharacterVector::create(),
    double sig_threshold = 0.05, double lfsr_null_threshold = 0.5) {
  // Extract feature IDs
  CharacterVector feature_ids = acat_matrix["feature_id"];
  int n_features = feature_ids.size();

  // Build bulk set for column filtering
  std::unordered_set<std::string> bulk_set;
  for (int i = 0; i < bulk_cts.size(); i++) {
    bulk_set.insert(as<std::string>(bulk_cts[i]));
  }

  // Get cell type columns (exclude feature_id and bulk columns)
  CharacterVector col_names = acat_matrix.names();
  std::vector<std::string> cell_types;
  for (int i = 0; i < col_names.size(); i++) {
    std::string col = as<std::string>(col_names[i]);
    if (col != "feature_id" && !bulk_set.count(col)) {
      cell_types.push_back(col);
    }
  }

  // Build L2 to L1 mapping if provided
  std::unordered_map<std::string, std::string> l2_l1_map =
      cascade::build_l2_to_l1_mapping(l2_to_l1_mapping);

  // Build feature->L1CT->lfsr map from long DataFrame
  std::unordered_map<std::string, std::unordered_map<std::string, double>>
      feature_lfsr_map;
  bool has_lfsr_lookup = lfsr_lookup.isNotNull();
  if (has_lfsr_lookup) {
    DataFrame lfsr_df = as<DataFrame>(lfsr_lookup);
    CharacterVector fids = lfsr_df["feature_id"];
    CharacterVector cts = lfsr_df["cell_type"];
    NumericVector vals = lfsr_df["min_lfsr"];
    for (int i = 0; i < fids.size(); i++) {
      double val = vals[i];
      if (NumericVector::is_na(val) || std::isnan(val)) continue;
      std::string fid = as<std::string>(fids[i]);
      std::string ct = as<std::string>(cts[i]);
      auto& ct_map = feature_lfsr_map[fid];
      auto it = ct_map.find(ct);
      if (it == ct_map.end() || val < it->second) {
        ct_map[ct] = val;
      }
    }
  }

  // Build generalized cell type sets
  cascade::CellTypeSets ct_sets(lineage_groups, subgroup_levels, bulk_cts,
                                other_cts);

  // Convert specificity categories to vector
  std::vector<std::string> categories =
      cascade::to_string_vector(specificity_categories);

  // Check for bulk columns in the data
  std::vector<std::string> bulk_col_names;
  std::vector<NumericVector> bulk_col_data;
  for (int i = 0; i < col_names.size(); i++) {
    std::string col = as<std::string>(col_names[i]);
    if (bulk_set.count(col)) {
      bulk_col_names.push_back(col);
      bulk_col_data.push_back(as<NumericVector>(acat_matrix[col]));
    }
  }

  // Pre-extract all cell type columns ONCE (avoid repeated DataFrame name
  // lookup)
  int n_ct = cell_types.size();
  std::vector<NumericVector> ct_columns(n_ct);
  for (int j = 0; j < n_ct; j++) {
    ct_columns[j] = as<NumericVector>(acat_matrix[cell_types[j]]);
  }

  // Prepare output vectors
  CharacterVector out_feature_ids(n_features);
  CharacterVector primary_categories(n_features);
  CharacterVector cell_type_specificity_pattern(n_features);
  CharacterVector significant_cts_str(n_features);
  CharacterVector tested_cts_str(n_features);
  IntegerVector n_significant_cts(n_features);
  IntegerVector n_tested_cts(n_features);

  // Process each feature
  for (int i = 0; i < n_features; i++) {
    std::string feature_id = as<std::string>(feature_ids[i]);
    out_feature_ids[i] = feature_id;

    // Collect significant and tested cell types
    std::vector<std::string> sig_cts_orig;
    std::vector<std::string> tested_cts;
    std::unordered_map<std::string, double> ct_qvalues;

    for (int j = 0; j < n_ct; j++) {
      double qval = ct_columns[j][i];
      if (!NumericVector::is_na(qval)) {
        tested_cts.push_back(cell_types[j]);
        ct_qvalues[cell_types[j]] = qval;
        if (qval < sig_threshold) {
          sig_cts_orig.push_back(cell_types[j]);
        }
      }
    }

    // Check bulk columns
    for (size_t bi = 0; bi < bulk_col_names.size(); bi++) {
      if (!NumericVector::is_na(bulk_col_data[bi][i])) {
        tested_cts.push_back(bulk_col_names[bi]);
        if (bulk_col_data[bi][i] < sig_threshold) {
          sig_cts_orig.push_back(bulk_col_names[bi]);
        }
      }
    }

    // Map to L1 and consolidate for categorization
    std::unordered_map<std::string, std::vector<std::string>> l1_to_l2_sigs;
    std::unordered_map<std::string, double> l1_min_qvals;
    std::vector<std::string> sig_cts_l1;

    for (const auto& ct : sig_cts_orig) {
      std::string l1_ct = cascade::map_to_l1(ct, l2_l1_map);

      // Track L2 to L1 mapping
      if (!l1_to_l2_sigs.count(l1_ct)) {
        l1_to_l2_sigs[l1_ct] = std::vector<std::string>();
        l1_min_qvals[l1_ct] = 1.0;
      }
      l1_to_l2_sigs[l1_ct].push_back(ct);

      // Update minimum q-value for L1
      if (ct_qvalues.count(ct) && ct_qvalues[ct] < l1_min_qvals[l1_ct]) {
        l1_min_qvals[l1_ct] = ct_qvalues[ct];
      }
    }

    // Get unique L1 significant cell types
    for (const auto& pair : l1_to_l2_sigs) {
      sig_cts_l1.push_back(pair.first);
    }

    // Build LFSR values for this feature
    std::unordered_map<std::string, double> lfsr_map;

    if (has_lfsr_lookup) {
      auto feat_it = feature_lfsr_map.find(feature_id);
      if (feat_it != feature_lfsr_map.end()) {
        // Build set of tested L1 CTs for this feature
        std::unordered_set<std::string> tested_l1_set;
        for (const auto& ct : tested_cts) {
          std::string l1_ct = cascade::map_to_l1(ct, l2_l1_map);
          tested_l1_set.insert(l1_ct);
        }
        // Filter to tested L1 CTs only
        for (const auto& pair : feat_it->second) {
          if (tested_l1_set.count(pair.first)) {
            lfsr_map[pair.first] = pair.second;
          }
        }
      }
    }

    // Categorize using generalized internal function
    std::string category = cascade::categorize_single_feature_internal(
        sig_cts_l1, lfsr_map, ct_sets, sig_threshold, lfsr_null_threshold,
        categories);

    primary_categories[i] = category;

    // Map to pattern number
    for (size_t j = 0; j < categories.size(); j++) {
      if (category == categories[j]) {
        cell_type_specificity_pattern[i] = std::to_string(j + 1);
        break;
      }
    }

    // Store results
    n_significant_cts[i] = sig_cts_orig.size();
    n_tested_cts[i] = tested_cts.size();

    // Join significant cell types (original, including L2)
    if (sig_cts_orig.empty()) {
      significant_cts_str[i] = "";
    } else {
      std::string sig_str = sig_cts_orig[0];
      for (size_t j = 1; j < sig_cts_orig.size(); j++) {
        sig_str += "," + sig_cts_orig[j];
      }
      significant_cts_str[i] = sig_str;
    }

    // Join tested cell types
    if (tested_cts.empty()) {
      tested_cts_str[i] = "";
    } else {
      std::string tested_str = tested_cts[0];
      for (size_t j = 1; j < tested_cts.size(); j++) {
        tested_str += "," + tested_cts[j];
      }
      tested_cts_str[i] = tested_str;
    }
  }

  // Build result DataFrame
  return DataFrame::create(
      Named("feature_id") = out_feature_ids,
      Named("primary_category") = primary_categories,
      Named("cell_type_specificity_pattern") = cell_type_specificity_pattern,
      Named("significant_cts") = significant_cts_str,
      Named("tested_cts") = tested_cts_str,
      Named("n_significant_cts") = n_significant_cts,
      Named("n_tested_cts") = n_tested_cts);
}

// [[Rcpp::export]]
CharacterVector categorize_variants_bulk_cpp(
    DataFrame variant_results, CharacterVector specificity_categories,
    Nullable<DataFrame> lfsr_lookup_results = R_NilValue,
    Nullable<List> l2_to_l1_mapping = R_NilValue,
    List lineage_groups = List::create(), List subgroup_levels = List::create(),
    CharacterVector bulk_cts = CharacterVector::create(),
    CharacterVector other_cts = CharacterVector::create(),
    double lfsr_sig_threshold = 0.05, double lfsr_null_threshold = 0.5) {
  // Extract columns from variant results
  CharacterVector variant_ids = variant_results["variant_id"];
  IntegerVector n_affected = variant_results["n_affected_cell_types"];
  CharacterVector affected_cts_str = variant_results["affected_cell_types"];

  int n_variants = variant_ids.size();
  CharacterVector specificities(n_variants);

  // Initialize all to last category (no significance)
  std::string nosig = as<std::string>(
      specificity_categories[specificity_categories.size() - 1]);
  for (int i = 0; i < n_variants; i++) {
    specificities[i] = nosig;
  }

  // Build L2 to L1 mapping
  std::unordered_map<std::string, std::string> l2_l1_map =
      cascade::build_l2_to_l1_mapping(l2_to_l1_mapping);

  // Build generalized cell type sets
  cascade::CellTypeSets ct_sets(lineage_groups, subgroup_levels, bulk_cts,
                                other_cts);

  // Convert specificity categories to vector
  std::vector<std::string> categories =
      cascade::to_string_vector(specificity_categories);

  // Process LFSR lookup results if provided
  std::unordered_map<std::string, std::unordered_map<std::string, double>>
      lfsr_data;
  if (lfsr_lookup_results.isNotNull()) {
    DataFrame lfsr_df = as<DataFrame>(lfsr_lookup_results);
    CharacterVector lfsr_var_ids = lfsr_df["variant_id"];
    CharacterVector lfsr_cell_types = lfsr_df["cell_type"];
    NumericVector lfsr_values = lfsr_df["min_lfsr"];

    for (int i = 0; i < lfsr_var_ids.size(); i++) {
      std::string var_id = as<std::string>(lfsr_var_ids[i]);
      std::string ct = as<std::string>(lfsr_cell_types[i]);
      double lfsr_val = lfsr_values[i];

      // Map to L1 if needed
      ct = cascade::map_to_l1(ct, l2_l1_map);

      // Store minimum LFSR for each variant-celltype pair
      if (!lfsr_data[var_id].count(ct) || lfsr_val < lfsr_data[var_id][ct]) {
        lfsr_data[var_id][ct] = lfsr_val;
      }
    }
  }

  // Process each variant
  for (int i = 0; i < n_variants; i++) {
    if (n_affected[i] == 0) continue;

    std::string var_id = as<std::string>(variant_ids[i]);

    // Parse affected cell types
    std::vector<std::string> sig_cts;
    std::string affected_str = as<std::string>(affected_cts_str[i]);

    if (!affected_str.empty() && affected_str != "NA") {
      // Split by comma and map to L1
      size_t pos = 0;
      while ((pos = affected_str.find(',')) != std::string::npos) {
        std::string ct = affected_str.substr(0, pos);
        // Map to L1 if needed
        if (l2_l1_map.count(ct)) {
          ct = l2_l1_map[ct];
        }
        // Add if not already present (unique L1)
        if (std::find(sig_cts.begin(), sig_cts.end(), ct) == sig_cts.end()) {
          sig_cts.push_back(ct);
        }
        affected_str.erase(0, pos + 1);
      }
      if (!affected_str.empty()) {
        std::string ct = affected_str;
        if (l2_l1_map.count(ct)) {
          ct = l2_l1_map[ct];
        }
        if (std::find(sig_cts.begin(), sig_cts.end(), ct) == sig_cts.end()) {
          sig_cts.push_back(ct);
        }
      }
    }

    // Get LFSR values for this variant
    std::unordered_map<std::string, double> lfsr_map;
    if (lfsr_data.count(var_id)) {
      lfsr_map = lfsr_data[var_id];
    }

    // Categorize using generalized internal function
    specificities[i] = cascade::categorize_single_feature_internal(
        sig_cts, lfsr_map, ct_sets, lfsr_sig_threshold, lfsr_null_threshold,
        categories);
  }

  return specificities;
}

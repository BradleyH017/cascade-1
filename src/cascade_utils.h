#ifndef CASCADE_UTILS_H
#define CASCADE_UTILS_H

#include <Rcpp.h>

#include <algorithm>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cascade {

// ============================================================================
// String Utilities
// ============================================================================

inline std::vector<std::string> split_string(const std::string& str,
                                             char delimiter) {
  std::vector<std::string> tokens;
  if (str.empty()) return tokens;

  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delimiter)) {
    if (!token.empty()) {
      tokens.push_back(token);
    }
  }
  return tokens;
}

inline std::string join_strings(const std::vector<std::string>& strings,
                                const std::string& delimiter) {
  if (strings.empty()) return "";

  std::stringstream ss;
  for (size_t i = 0; i < strings.size(); ++i) {
    if (i > 0) ss << delimiter;
    ss << strings[i];
  }
  return ss.str();
}

// ============================================================================
// Cell Type Set Utilities
// ============================================================================

// Generalized structure to hold cell type sets for categorization.
// Supports N lineages, M subgroup levels, and configurable bulk/other sets.
struct CellTypeSets {
  // Top-level lineage groups (N groups, N >= 2)
  std::vector<std::unordered_set<std::string>> lineages;

  // Subgroup levels ordered broadest → narrowest.
  // subgroup_levels[0] = broadest level (vector of groups)
  // subgroup_levels[1] = next narrower level, etc.
  std::vector<std::vector<std::unordered_set<std::string>>> subgroup_levels;

  // Bulk cell types (e.g., PBMC) — bulk-only signal → underpowered
  std::unordered_set<std::string> bulk;

  // "Other" cell types — excluded from lineage counting but participate
  // in single-cell-type detection
  std::unordered_set<std::string> other;

  // Constructor from R Lists
  CellTypeSets(const Rcpp::List& lineage_list,
               const Rcpp::List& subgroup_levels_list,
               const Rcpp::CharacterVector& bulk_cts,
               const Rcpp::CharacterVector& other_cts) {
    // Build lineage sets
    for (int i = 0; i < lineage_list.size(); i++) {
      Rcpp::CharacterVector cts =
          Rcpp::as<Rcpp::CharacterVector>(lineage_list[i]);
      std::unordered_set<std::string> s;
      for (int j = 0; j < cts.size(); j++) {
        s.insert(Rcpp::as<std::string>(cts[j]));
      }
      lineages.push_back(std::move(s));
    }
    // Build subgroup levels
    for (int i = 0; i < subgroup_levels_list.size(); i++) {
      Rcpp::List level = Rcpp::as<Rcpp::List>(subgroup_levels_list[i]);
      std::vector<std::unordered_set<std::string>> level_groups;
      for (int j = 0; j < level.size(); j++) {
        Rcpp::CharacterVector cts = Rcpp::as<Rcpp::CharacterVector>(level[j]);
        std::unordered_set<std::string> s;
        for (int k = 0; k < cts.size(); k++) {
          s.insert(Rcpp::as<std::string>(cts[k]));
        }
        level_groups.push_back(std::move(s));
      }
      subgroup_levels.push_back(std::move(level_groups));
    }
    // Build bulk and other sets
    for (int i = 0; i < bulk_cts.size(); i++) {
      bulk.insert(Rcpp::as<std::string>(bulk_cts[i]));
    }
    for (int i = 0; i < other_cts.size(); i++) {
      other.insert(Rcpp::as<std::string>(other_cts[i]));
    }
  }
};

// Build L2 to L1 mapping from R List
inline std::unordered_map<std::string, std::string> build_l2_to_l1_mapping(
    const Rcpp::Nullable<Rcpp::List>& l2_to_l1_mapping) {
  std::unordered_map<std::string, std::string> l2_l1_map;

  if (l2_to_l1_mapping.isNotNull()) {
    Rcpp::List mapping = Rcpp::as<Rcpp::List>(l2_to_l1_mapping);
    Rcpp::CharacterVector l2_names = mapping.names();
    for (int i = 0; i < l2_names.size(); i++) {
      std::string l2_ct = Rcpp::as<std::string>(l2_names[i]);
      std::string l1_ct = Rcpp::as<std::string>(mapping[i]);
      l2_l1_map[l2_ct] = l1_ct;
    }
  }

  return l2_l1_map;
}

// Build L2 to L1 mapping from parallel vectors
inline std::unordered_map<std::string, std::string> build_l2_to_l1_mapping(
    const Rcpp::CharacterVector& l2_keys,
    const Rcpp::CharacterVector& l1_values) {
  std::unordered_map<std::string, std::string> l2_l1_map;
  for (int i = 0; i < l2_keys.size(); i++) {
    l2_l1_map[Rcpp::as<std::string>(l2_keys[i])] =
        Rcpp::as<std::string>(l1_values[i]);
  }
  return l2_l1_map;
}

// Map cell type to L1 level
inline std::string map_to_l1(
    const std::string& cell_type,
    const std::unordered_map<std::string, std::string>& l2_l1_map) {
  auto it = l2_l1_map.find(cell_type);
  return (it != l2_l1_map.end()) ? it->second : cell_type;
}

// Convert CharacterVector to vector of strings
inline std::vector<std::string> to_string_vector(
    const Rcpp::CharacterVector& char_vec) {
  std::vector<std::string> result;
  result.reserve(char_vec.size());
  for (int i = 0; i < char_vec.size(); i++) {
    result.push_back(Rcpp::as<std::string>(char_vec[i]));
  }
  return result;
}

// ============================================================================
// Cell Type Specificity Categorization
// ============================================================================

// LFSR gray zone check with configurable scope.
// exclude_set = nullptr → check ALL non-sig non-bulk types (single-cell-type
// behavior)
//   NOTE: "other" types ARE checked in this case
// exclude_set = &group → skip types in that group AND "other" types
//   (subgroup/lineage behavior: "other" is outside all lineages)
inline bool has_lfsr_gray_zone(
    const std::unordered_map<std::string, double>& lfsr_map,
    const std::vector<std::string>& sig_cts, const CellTypeSets& ct_sets,
    const std::unordered_set<std::string>* exclude_set,
    double lfsr_sig_threshold, double lfsr_null_threshold) {
  if (lfsr_map.empty()) return false;
  std::unordered_set<std::string> sig_set(sig_cts.begin(), sig_cts.end());
  for (const auto& pair : lfsr_map) {
    const std::string& ct = pair.first;
    double val = pair.second;
    if (sig_set.count(ct)) continue;
    if (ct_sets.bulk.count(ct)) continue;
    if (exclude_set) {
      // Subgroup/lineage mode: also skip "other" types
      if (ct_sets.other.count(ct)) continue;
      if (exclude_set->count(ct)) continue;
    }
    // Single-cell-type mode (exclude_set=nullptr): check all including "other"
    if (val >= lfsr_sig_threshold && val < lfsr_null_threshold) return true;
  }
  return false;
}

// Internal function for categorizing cell type specificity.
// Generalized to N lineages + M subgroup levels.
// Categories layout: [cross-lineage, underpowered,
//   lineage-specific, subgroup_1-specific, ..., subgroup_M-specific,
//   single-cell-type, no-significance]
std::string categorize_single_feature_internal(
    const std::vector<std::string>& significant_cts,
    const std::unordered_map<std::string, double>& lfsr_map,
    const CellTypeSets& ct_sets, double lfsr_sig_threshold,
    double lfsr_null_threshold, const std::vector<std::string>& categories);

}  // namespace cascade

#endif  // CASCADE_UTILS_H

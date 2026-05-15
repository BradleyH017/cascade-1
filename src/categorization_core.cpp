#include "cascade_utils.h"

namespace cascade {

// Generalized categorization: N lineages + M subgroup levels.
// Categories layout: [cross-lineage, underpowered,
//   lineage-specific, subgroup_1-specific, ..., subgroup_M-specific,
//   single-cell-type, no-significance]
std::string categorize_single_feature_internal(
    const std::vector<std::string>& significant_cts,
    const std::unordered_map<std::string, double>& lfsr_map,
    const CellTypeSets& ct_sets, double lfsr_sig_threshold,
    double lfsr_null_threshold, const std::vector<std::string>& categories) {
  const int idx_cross = 0;
  const int idx_underpow = 1;
  const int idx_lineage = 2;
  // Subgroup categories: indices 3 .. (3 + n_subgroup_levels - 1)
  const int idx_single = static_cast<int>(categories.size()) - 2;
  const int idx_nosig = static_cast<int>(categories.size()) - 1;

  // Remove bulk from sig types. "other" types stay in sig_core —
  // they participate in single-cell-type detection.
  std::vector<std::string> sig_core;
  bool has_bulk = false;
  for (const auto& ct : significant_cts) {
    if (ct_sets.bulk.count(ct)) {
      has_bulk = true;
      continue;
    }
    sig_core.push_back(ct);
  }

  // Bulk-only → underpowered
  if (sig_core.empty() && has_bulk) return categories[idx_underpow];
  // Nothing significant
  if (sig_core.empty()) return categories[idx_nosig];

  // Count lineages with signal.
  // "other" types are not in any lineage set, so naturally skipped.
  int n_sig_lineages = 0;
  std::vector<bool> lineage_has_sig(ct_sets.lineages.size(), false);
  for (size_t li = 0; li < ct_sets.lineages.size(); li++) {
    for (const auto& ct : sig_core) {
      if (ct_sets.lineages[li].count(ct)) {
        lineage_has_sig[li] = true;
        n_sig_lineages++;
        break;
      }
    }
  }

  // Cross-lineage: NO LFSR check (broadest category; cannot be demoted)
  if (n_sig_lineages >= 2) {
    return categories[idx_cross];
  }

  // Single cell type: LFSR checks ALL non-sig, non-bulk tested types
  if (sig_core.size() == 1) {
    if (has_lfsr_gray_zone(lfsr_map, sig_core, ct_sets,
                           nullptr,  // check everything including "other"
                           lfsr_sig_threshold, lfsr_null_threshold))
      return categories[idx_underpow];
    return categories[idx_single];
  }

  // Check subgroup levels from NARROWEST to BROADEST
  for (int lev = static_cast<int>(ct_sets.subgroup_levels.size()) - 1; lev >= 0;
       lev--) {
    for (const auto& group : ct_sets.subgroup_levels[lev]) {
      bool all_in = true;
      for (const auto& ct : sig_core) {
        if (!group.count(ct)) {
          all_in = false;
          break;
        }
      }
      if (all_in) {
        int cat_idx = 3 + lev;
        // LFSR: check sibling groups + other lineages (NOT within-group types)
        if (has_lfsr_gray_zone(lfsr_map, sig_core, ct_sets,
                               &group,  // exclude same-group types and "other"
                               lfsr_sig_threshold, lfsr_null_threshold))
          return categories[idx_underpow];
        return categories[cat_idx];
      }
    }
  }

  // Lineage-specific: LFSR checks OTHER lineages only
  if (n_sig_lineages == 1) {
    // Build set of types in the significant lineage
    const std::unordered_set<std::string>* sig_lineage = nullptr;
    for (size_t li = 0; li < ct_sets.lineages.size(); li++) {
      if (lineage_has_sig[li]) {
        sig_lineage = &ct_sets.lineages[li];
        break;
      }
    }
    if (sig_lineage) {
      if (has_lfsr_gray_zone(
              lfsr_map, sig_core, ct_sets,
              sig_lineage,  // exclude same-lineage types and "other"
              lfsr_sig_threshold, lfsr_null_threshold))
        return categories[idx_underpow];
    }
    return categories[idx_lineage];
  }

  return categories[idx_nosig];
}

}  // namespace cascade

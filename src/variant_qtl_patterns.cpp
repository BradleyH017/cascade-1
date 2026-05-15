#include <Rcpp.h>

#include <algorithm>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace Rcpp;

// Helper function to extract character vector from 2D list array
CharacterVector getArrayElement(RObject arr, int row, int col) {
  if (arr.isNULL()) return CharacterVector();

  try {
    // R stores matrices column-wise, so the linear index is: row + col * nrow
    IntegerVector dims = arr.attr("dim");
    if (dims.size() != 2) return CharacterVector();

    int nrow = dims[0];
    int idx = row + col * nrow;

    // Access the underlying list
    List arr_list = as<List>(arr);
    if (idx >= arr_list.size()) return CharacterVector();

    SEXP elem = arr_list[idx];
    if (Rf_isNull(elem)) return CharacterVector();
    if (TYPEOF(elem) != STRSXP) return CharacterVector();

    return as<CharacterVector>(elem);
  } catch (...) {
    return CharacterVector();
  }
}

// [[Rcpp::export]]
List detect_variant_qtl_patterns(
    CharacterVector variant_ids, List variant_matrix_data,
    DataFrame peak_gene_links, CharacterVector cell_types, int ct_idx,
    CharacterVector qtl_mechanism_categories = CharacterVector::create(
        "Local Cascade", "Positional Cascade", "Distal Cascade",
        "caQTL + eQTL (No Link)", "Only caQTL (With Link)",
        "Only caQTL (No Link)", "Only eQTL", "No molQTL"),
    IntegerVector pattern_to_mechanism = IntegerVector::create(
        1, 2, 2, 2, 3, 4, 4, 4, 4, 4,  // Patterns 1-10
        4, 4, 5, 5, 5, 5, 6, 6, 6, 7,  // Patterns 11-20
        7, 7, 8, 8, 8                  // Patterns 21-25
        ),
    List caqtl_status_values =
        List::create(Named("OVERLAPPING") = "For overlapping peak",
                     Named("NON_OVERLAPPING") = "For non-overlapping peak",
                     Named("NO_CAQTL") = "No detected caQTL"),
    List eqtl_status_values = List::create(
        Named("LINKED") = "For linked gene",
        Named("NON_LINKED") = "For non-linked gene", Named("EQTL") = "eQTL",
        Named("NO_EQTL") = "No detected eQTL")) {
  int n_variants = variant_ids.size();

  // Extract matrices and lists
  LogicalMatrix caqtl_matrix = variant_matrix_data["caqtl_matrix"];
  LogicalMatrix eqtl_matrix = variant_matrix_data["eqtl_matrix"];
  List peak_overlaps = variant_matrix_data["peak_overlaps"];

  // Extract the peak and gene lists
  RObject caqtl_peaks_list = variant_matrix_data["caqtl_peaks_list"];
  RObject eqtl_genes_list = variant_matrix_data["eqtl_genes_list"];

  // Get cell-type specific masks
  LogicalVector caqtl_mask = caqtl_matrix(_, ct_idx);
  LogicalVector eqtl_mask = eqtl_matrix(_, ct_idx);

  // Initialize result vectors
  IntegerVector qtl_patterns(n_variants);
  CharacterVector qtl_categories(n_variants);
  // New columns to replace cascade_gene and affected_peaks
  CharacterVector associated_genes_vec(n_variants, NA_STRING);
  CharacterVector associated_peaks_vec(n_variants, NA_STRING);
  CharacterVector cascade_peak_genes_vec(n_variants, NA_STRING);

  // Initialize new column vectors
  LogicalVector peak_overlap(n_variants, false);
  CharacterVector caqtl_status(n_variants);
  CharacterVector peak_gene_link_status(
      n_variants);  // Changed to CharacterVector for categorical values
  CharacterVector eqtl_status(n_variants);
  CharacterVector link_mechanism_vec(n_variants, NA_STRING);

  // Create peak-gene lookup and mechanism lookup
  std::unordered_map<std::string, std::vector<std::string>> peak_to_genes;
  // Key: "peak_id::gene_id" -> mechanism string
  std::unordered_map<std::string, std::string> peak_gene_to_mechanism;
  CharacterVector peak_ids = peak_gene_links["peak_id"];
  CharacterVector gene_ids = peak_gene_links["gene_id"];

  // Check if mechanism column exists
  bool has_mechanism_col = peak_gene_links.containsElementNamed("mechanism");
  CharacterVector mechanism_col;
  if (has_mechanism_col) {
    mechanism_col = peak_gene_links["mechanism"];
  }

  for (int i = 0; i < peak_ids.size(); i++) {
    std::string pid = as<std::string>(peak_ids[i]);
    std::string gid = as<std::string>(gene_ids[i]);
    peak_to_genes[pid].push_back(gid);
    if (has_mechanism_col) {
      peak_gene_to_mechanism[pid + "::" + gid] =
          as<std::string>(mechanism_col[i]);
    }
  }

  // Process each variant with complete pattern detection
  for (int i = 0; i < n_variants; i++) {
    bool has_caqtl = caqtl_mask[i];
    bool has_eqtl = eqtl_mask[i];

    // Get overlapping peaks for this variant
    bool has_overlap = false;
    CharacterVector overlaps;

    // Get the current cell type
    std::string current_cell_type = as<std::string>(cell_types[ct_idx]);

    // Safely extract overlaps from the list
    if (i < peak_overlaps.size()) {
      SEXP overlap_elem = peak_overlaps[i];
      if (!Rf_isNull(overlap_elem)) {
        if (TYPEOF(overlap_elem) == STRSXP) {
          // Error: peak overlaps must be a list with cell type information
          std::string variant_name = as<std::string>(variant_ids[i]);
          Rcpp::stop("Peak overlap data for variant " + variant_name +
                     " must be a list with cell_type information; got a "
                     "character vector.");
        } else if (TYPEOF(overlap_elem) == VECSXP) {
          // It's a list of peak info - extract peak_ids only for the current
          // cell type
          List overlap_list = as<List>(overlap_elem);
          if (overlap_list.size() > 0) {
            // Extract peak_ids from the list, filtering by cell type
            CharacterVector peak_ids;
            for (int j = 0; j < overlap_list.size(); j++) {
              SEXP item = overlap_list[j];
              if (TYPEOF(item) == VECSXP) {
                List peak_info = as<List>(item);
                // Check if this peak is for the current cell type
                if (peak_info.containsElementNamed("cell_type") &&
                    peak_info.containsElementNamed("peak_id")) {
                  CharacterVector ct_vec = peak_info["cell_type"];
                  if (ct_vec.size() > 0 &&
                      as<std::string>(ct_vec[0]) == current_cell_type) {
                    CharacterVector pid = peak_info["peak_id"];
                    if (pid.size() > 0) {
                      peak_ids.push_back(as<std::string>(pid[0]));
                    }
                  }
                } else if (peak_info.containsElementNamed("peak_id")) {
                  // Error: peak overlaps must have cell type information
                  std::string variant_name = as<std::string>(variant_ids[i]);
                  Rcpp::stop("Peak overlap data for variant " + variant_name +
                             " is missing cell_type information. All peak "
                             "overlaps must specify the cell type.");
                }
              }
            }
            overlaps = peak_ids;
            has_overlap = overlaps.size() > 0;
          }
        }
      }
    }

    // Initialize pattern detection variables
    bool has_caqtl_overlapping = false;
    bool has_caqtl_nonoverlapping = false;
    bool has_link_overlapping = false;
    bool has_link_nonoverlapping = false;
    bool eqtl_for_linked_overlapping = false;
    bool eqtl_for_linked_nonoverlapping = false;

    std::vector<std::string> all_linked_genes;
    std::vector<std::string> overlapping_affected_peaks;
    std::vector<std::string> nonoverlapping_affected_peaks;

    // Check if overlapping peaks have gene links (needed for Pattern 11)
    if (has_overlap) {
      for (int j = 0; j < overlaps.size(); j++) {
        std::string peak = as<std::string>(overlaps[j]);
        if (peak_to_genes.find(peak) != peak_to_genes.end()) {
          has_link_overlapping = true;
          // Don't add to all_linked_genes here - will be done in caQTL section
          // if needed
          break;
        }
      }
    }

    // Determine caQTL status for overlapping and non-overlapping peaks
    if (has_caqtl) {
      // Get actual caQTL peaks for this variant
      CharacterVector caqtl_peaks;

      // Extract from the 2D array-like structure
      if (!caqtl_peaks_list.isNULL()) {
        try {
          // R stores matrices column-wise, so index = row + col * nrow
          IntegerVector dims = caqtl_peaks_list.attr("dim");
          if (dims.size() == 2) {
            int nrow = dims[0];
            int idx = i + ct_idx * nrow;

            List peaks_as_list = as<List>(caqtl_peaks_list);
            if (idx < peaks_as_list.size()) {
              SEXP elem = peaks_as_list[idx];
              if (!Rf_isNull(elem) && TYPEOF(elem) == STRSXP) {
                caqtl_peaks = as<CharacterVector>(elem);
              }
            }
          }
        } catch (...) {
          // If extraction fails, leave caqtl_peaks empty
        }
      }

      // Separate caQTL peaks into overlapping and non-overlapping
      if (has_overlap && caqtl_peaks.size() > 0) {
        // Check which caQTL peaks overlap with the variant
        std::unordered_set<std::string> overlap_set;
        for (int j = 0; j < overlaps.size(); j++) {
          overlap_set.insert(as<std::string>(overlaps[j]));
        }

        for (int j = 0; j < caqtl_peaks.size(); j++) {
          std::string peak = as<std::string>(caqtl_peaks[j]);
          if (overlap_set.find(peak) != overlap_set.end()) {
            // This caQTL peak overlaps with the variant
            has_caqtl_overlapping = true;
            overlapping_affected_peaks.push_back(peak);

            if (peak_to_genes.find(peak) != peak_to_genes.end()) {
              has_link_overlapping = true;
              for (const auto& gene : peak_to_genes[peak]) {
                all_linked_genes.push_back(gene);
              }
            }
          } else {
            // This caQTL peak doesn't overlap with the variant
            has_caqtl_nonoverlapping = true;
            nonoverlapping_affected_peaks.push_back(peak);

            if (peak_to_genes.find(peak) != peak_to_genes.end()) {
              has_link_nonoverlapping = true;
              for (const auto& gene : peak_to_genes[peak]) {
                all_linked_genes.push_back(gene);
              }
            }
          }
        }

        // If we still have no peaks but has_caqtl is true, use fallback
        if (overlapping_affected_peaks.empty() &&
            nonoverlapping_affected_peaks.empty() && caqtl_peaks.size() == 0) {
          has_caqtl_nonoverlapping = true;
          nonoverlapping_affected_peaks.push_back("distal_peak");
        }
      } else {
        // No overlap but has caQTL - all peaks are non-overlapping
        has_caqtl_nonoverlapping = true;

        // Process all caQTL peaks as non-overlapping
        for (int j = 0; j < caqtl_peaks.size(); j++) {
          std::string peak = as<std::string>(caqtl_peaks[j]);
          nonoverlapping_affected_peaks.push_back(peak);

          // Check if this peak has gene links
          if (peak_to_genes.find(peak) != peak_to_genes.end()) {
            has_link_nonoverlapping = true;
            for (const auto& gene : peak_to_genes[peak]) {
              all_linked_genes.push_back(gene);
            }
          }
        }

        // If no peaks found, use "distal_peak" as fallback
        if (nonoverlapping_affected_peaks.empty()) {
          nonoverlapping_affected_peaks.push_back("distal_peak");
        }
      }
    }

    // Extract eQTL genes for this variant
    std::vector<std::string> eqtl_genes_for_variant;
    if (has_eqtl && !eqtl_genes_list.isNULL()) {
      try {
        int eqtl_idx = i + ct_idx * n_variants;
        SEXP eqtl_elem = VECTOR_ELT(eqtl_genes_list, eqtl_idx);
        if (!Rf_isNull(eqtl_elem) && TYPEOF(eqtl_elem) == STRSXP) {
          CharacterVector genes = as<CharacterVector>(eqtl_elem);
          for (int j = 0; j < genes.size(); j++) {
            eqtl_genes_for_variant.push_back(as<std::string>(genes[j]));
          }
        }
      } catch (...) {
        // If extraction fails, continue with empty eqtl_genes
      }
    }

    // Check if overlapping peaks (independent of caQTL status) link to eQTL
    // genes This is needed for patterns 4, 5, and others
    bool has_overlap_peak_link = false;
    bool eqtl_for_overlap_peak_link = false;

    if (has_overlap && !overlaps.isNULL() && overlaps.size() > 0) {
      // Check if any overlapping peak has a gene link
      for (int j = 0; j < overlaps.size(); j++) {
        std::string peak = as<std::string>(overlaps[j]);
        if (peak_to_genes.find(peak) != peak_to_genes.end() &&
            !peak_to_genes[peak].empty()) {
          has_overlap_peak_link = true;

          // If we have eQTL genes, check if they match
          if (has_eqtl && !eqtl_genes_for_variant.empty()) {
            std::set<std::string> eqtl_set(eqtl_genes_for_variant.begin(),
                                           eqtl_genes_for_variant.end());
            for (const auto& gene : peak_to_genes[peak]) {
              if (eqtl_set.find(gene) != eqtl_set.end()) {
                eqtl_for_overlap_peak_link = true;
                all_linked_genes.push_back(gene);
                break;
              }
            }
          }
        }
        if (eqtl_for_overlap_peak_link) break;
      }
    }

    // Pattern-3 path: variant overlaps a peak with a peak-gene link, has an
    // eQTL for the linked gene, but no caQTL. Mark the eQTL as linked via the
    // overlapping peak so downstream classification reaches Positional Cascade.
    if (!has_caqtl && has_overlap && has_overlap_peak_link && has_eqtl &&
        eqtl_for_overlap_peak_link) {
      eqtl_for_linked_overlapping = true;
    }

    // Check if eQTL affects linked genes (for caQTL patterns)
    if (has_caqtl && has_eqtl && !all_linked_genes.empty() &&
        !eqtl_genes_for_variant.empty()) {
      // Check actual overlap between eQTL genes and linked genes
      std::set<std::string> linked_set(all_linked_genes.begin(),
                                       all_linked_genes.end());
      std::set<std::string> eqtl_set(eqtl_genes_for_variant.begin(),
                                     eqtl_genes_for_variant.end());

      // Find intersection
      std::vector<std::string> intersection;
      std::set_intersection(linked_set.begin(), linked_set.end(),
                            eqtl_set.begin(), eqtl_set.end(),
                            std::back_inserter(intersection));

      if (!intersection.empty()) {
        // There are genes that are both linked and eQTL targets
        if (has_link_overlapping) {
          // Check if any overlapping peaks link to eQTL genes
          for (const auto& peak : overlapping_affected_peaks) {
            if (peak_to_genes.find(peak) != peak_to_genes.end()) {
              for (const auto& gene : peak_to_genes[peak]) {
                if (eqtl_set.find(gene) != eqtl_set.end()) {
                  eqtl_for_linked_overlapping = true;
                  break;
                }
              }
            }
            if (eqtl_for_linked_overlapping) break;
          }
        }
        if (has_link_nonoverlapping) {
          // Check if any non-overlapping peaks link to eQTL genes
          for (const auto& peak : nonoverlapping_affected_peaks) {
            if (peak_to_genes.find(peak) != peak_to_genes.end()) {
              for (const auto& gene : peak_to_genes[peak]) {
                if (eqtl_set.find(gene) != eqtl_set.end()) {
                  eqtl_for_linked_nonoverlapping = true;
                  break;
                }
              }
            }
            if (eqtl_for_linked_nonoverlapping) break;
          }
        }
      }
    }

    // Apply 25-pattern MECE logic
    int pattern = 25;  // Default: No molQTL (no overlap)

    if (has_overlap) {
      if (has_caqtl_overlapping) {
        // Variant overlaps peak and is caQTL for that overlapping peak
        if (has_link_overlapping) {
          if (has_eqtl && eqtl_for_linked_overlapping) {
            pattern = 1;  // Local Cascade
          } else if (has_eqtl && !eqtl_for_linked_overlapping) {
            pattern = 6;  // caQTL + eQTL (overlap link, discordant)
          } else {
            pattern = 13;  // Only caQTL (overlap link)
          }
        } else {
          if (has_eqtl) {
            pattern = 10;  // caQTL + eQTL (no link)
          } else {
            pattern = 17;  // Only caQTL (overlap, no link)
          }
        }
      } else if (has_caqtl_nonoverlapping) {
        // Variant overlaps peak but is caQTL for a different (non-overlapping)
        // peak Check overlapping-peak link first (higher priority than
        // non-overlap link)
        if (has_overlap_peak_link) {
          if (has_eqtl && eqtl_for_overlap_peak_link) {
            pattern =
                2;  // Positional Cascade (overlap link, non-overlap caQTL)
          } else if (has_eqtl && !eqtl_for_overlap_peak_link) {
            pattern = 7;  // caQTL + eQTL (overlap link, non-overlap caQTL,
                          // discordant)
          } else {
            pattern = 14;  // Only caQTL (overlap link, non-overlap caQTL)
          }
        } else if (has_link_nonoverlapping) {
          // Link via non-overlapping caQTL peak
          if (has_eqtl && eqtl_for_linked_nonoverlapping) {
            pattern = 4;  // Positional Cascade (non-overlap link)
          } else if (has_eqtl && !eqtl_for_linked_nonoverlapping) {
            pattern = 8;  // caQTL + eQTL (non-overlap link, discordant)
          } else {
            pattern = 15;  // Only caQTL (non-overlap link)
          }
        } else {
          // No links from either overlapping or non-overlapping peaks
          if (has_eqtl) {
            pattern = 11;  // caQTL + eQTL (non-overlap caQTL, no link)
          } else {
            pattern = 18;  // Only caQTL (non-overlap caQTL, no link)
          }
        }
      } else {
        // Overlap, but no caQTL
        if (has_overlap_peak_link) {
          if (has_eqtl && eqtl_for_overlap_peak_link) {
            pattern = 3;  // Positional Cascade (overlap link)
          } else if (has_eqtl && !eqtl_for_overlap_peak_link) {
            pattern = 20;  // Only eQTL (overlap link)
          } else {
            pattern = 23;  // No molQTL (overlap link)
          }
        } else {
          if (has_eqtl) {
            pattern = 21;  // Only eQTL (no link)
          } else {
            pattern = 24;  // No molQTL (no link)
          }
        }
      }
    } else {
      // No overlap cases
      if (has_caqtl_nonoverlapping) {
        if (has_link_nonoverlapping) {
          if (has_eqtl && eqtl_for_linked_nonoverlapping) {
            pattern = 5;  // Distal Cascade
          } else if (has_eqtl && !eqtl_for_linked_nonoverlapping) {
            pattern = 9;  // caQTL + eQTL (no overlap, discordant)
          } else {
            pattern = 16;  // Only caQTL (no overlap, link)
          }
        } else {
          if (has_eqtl) {
            pattern = 12;  // caQTL + eQTL (no overlap, no link)
          } else {
            pattern = 19;  // Only caQTL (no overlap, no link)
          }
        }
      } else {
        // No caQTL and no overlap
        if (has_eqtl) {
          pattern = 22;  // Only eQTL (no overlap)
        } else {
          pattern = 25;  // No molQTL (no overlap)
        }
      }
    }

    qtl_patterns[i] = pattern;

    // Set new column values based on pattern and data
    peak_overlap[i] = has_overlap;

    // Set caQTL status using passed constants
    if (has_caqtl_overlapping) {
      caqtl_status[i] = as<std::string>(caqtl_status_values["OVERLAPPING"]);
    } else if (has_caqtl_nonoverlapping) {
      caqtl_status[i] = as<std::string>(caqtl_status_values["NON_OVERLAPPING"]);
    } else {
      caqtl_status[i] = as<std::string>(caqtl_status_values["NO_CAQTL"]);
    }

    // Set peak-gene link status from pattern number.
    // "For caQTL peak": link is from the caQTL-linked peak (same peak as caQTL)
    //   patterns: 1, 4, 5, 6, 8, 9, 13, 15, 16
    // "For overlapping peak": link is from the variant's overlapping peak, but
    // caQTL
    //   is either from a different peak or absent
    //   patterns: 2, 3, 7, 14, 20, 23
    // "No link": all others
    if (pattern == 1 || pattern == 4 || pattern == 5 || pattern == 6 ||
        pattern == 8 || pattern == 9 || pattern == 13 || pattern == 15 ||
        pattern == 16) {
      peak_gene_link_status[i] = "For caQTL peak";
    } else if (pattern == 2 || pattern == 3 || pattern == 7 || pattern == 14 ||
               pattern == 20 || pattern == 23) {
      peak_gene_link_status[i] = "For overlapping peak";
    } else {
      peak_gene_link_status[i] = "No link";
    }

    // Set eQTL status from pattern number.
    if (has_eqtl) {
      // LINKED — cascade patterns where eQTL matches the linked gene
      //   patterns: 1, 2, 3, 4, 5 (all cascade)
      if (pattern == 1 || pattern == 2 || pattern == 3 || pattern == 4 ||
          pattern == 5) {
        eqtl_status[i] = as<std::string>(eqtl_status_values["LINKED"]);
      } else if (pattern == 6 || pattern == 7 || pattern == 8 || pattern == 9 ||
                 pattern == 20) {
        // NON_LINKED — discordant: link exists but eQTL targets different gene
        //   patterns: 6, 7, 8, 9, 20
        eqtl_status[i] = as<std::string>(eqtl_status_values["NON_LINKED"]);
      } else {
        // EQTL — eQTL exists but no peak-gene link
        //   patterns: 10, 11, 12, 21, 22
        eqtl_status[i] = as<std::string>(eqtl_status_values["EQTL"]);
      }
    } else {
      eqtl_status[i] = as<std::string>(eqtl_status_values["NO_EQTL"]);
    }

    // Collect associated genes (all eQTL genes)
    std::unordered_set<std::string> all_associated_genes;
    if (has_eqtl && !eqtl_genes_for_variant.empty()) {
      all_associated_genes.insert(eqtl_genes_for_variant.begin(),
                                  eqtl_genes_for_variant.end());
    }

    // Collect cascade peak->gene pairs for cascade patterns (NEW: 1 Local, 2-4
    // Positional, 5 Distal)
    std::vector<std::string> cascade_peak_gene_pairs;
    if ((pattern == 1 || pattern == 2 || pattern == 3 || pattern == 4 ||
         pattern == 5) &&
        !all_linked_genes.empty() && has_eqtl &&
        !eqtl_genes_for_variant.empty()) {
      // Find genes that are both linked and eQTL targets
      std::set<std::string> linked_set(all_linked_genes.begin(),
                                       all_linked_genes.end());
      std::set<std::string> eqtl_set(eqtl_genes_for_variant.begin(),
                                     eqtl_genes_for_variant.end());

      std::vector<std::string> cascade_genes;
      std::set_intersection(linked_set.begin(), linked_set.end(),
                            eqtl_set.begin(), eqtl_set.end(),
                            std::back_inserter(cascade_genes));

      // Now find which peaks link to these cascade genes
      if (!cascade_genes.empty()) {
        // For Positional Cascade patterns where the gene link is via the
        // overlapping peak (P2, P3), only check overlapping peaks. The caQTL
        // peak (if any) is a different non-overlapping peak whose link is not
        // the cascade link.
        if ((pattern == 2 || pattern == 3) && has_overlap &&
            !overlaps.isNULL()) {
          for (int j = 0; j < overlaps.size(); j++) {
            std::string peak = as<std::string>(overlaps[j]);
            if (peak_to_genes.find(peak) != peak_to_genes.end()) {
              for (const auto& gene : peak_to_genes[peak]) {
                if (std::find(cascade_genes.begin(), cascade_genes.end(),
                              gene) != cascade_genes.end()) {
                  cascade_peak_gene_pairs.push_back(peak + "->" + gene);
                }
              }
            }
          }
        } else {
          // For Local Cascade (P1), Positional Cascade non-overlap link (P4),
          // and Distal Cascade (P5): check both overlapping and non-overlapping
          // caQTL peaks.
          for (const auto& peak : overlapping_affected_peaks) {
            if (peak_to_genes.find(peak) != peak_to_genes.end()) {
              for (const auto& gene : peak_to_genes[peak]) {
                if (std::find(cascade_genes.begin(), cascade_genes.end(),
                              gene) != cascade_genes.end()) {
                  cascade_peak_gene_pairs.push_back(peak + "->" + gene);
                }
              }
            }
          }
          for (const auto& peak : nonoverlapping_affected_peaks) {
            if (peak_to_genes.find(peak) != peak_to_genes.end()) {
              for (const auto& gene : peak_to_genes[peak]) {
                if (std::find(cascade_genes.begin(), cascade_genes.end(),
                              gene) != cascade_genes.end()) {
                  cascade_peak_gene_pairs.push_back(peak + "->" + gene);
                }
              }
            }
          }
        }
      }
    }

    // Combine all affected peaks
    std::vector<std::string> all_affected_peaks;
    all_affected_peaks.insert(all_affected_peaks.end(),
                              overlapping_affected_peaks.begin(),
                              overlapping_affected_peaks.end());
    all_affected_peaks.insert(all_affected_peaks.end(),
                              nonoverlapping_affected_peaks.begin(),
                              nonoverlapping_affected_peaks.end());

    // Set associated_genes
    if (!all_associated_genes.empty()) {
      std::string genes_str;
      size_t j = 0;
      for (const auto& gene : all_associated_genes) {
        if (j > 0) genes_str += ",";
        genes_str += gene;
        j++;
      }
      associated_genes_vec[i] = genes_str;
    }

    // Set associated_peaks
    if (!all_affected_peaks.empty()) {
      // Remove duplicates by using a set
      std::unordered_set<std::string> unique_peaks(all_affected_peaks.begin(),
                                                   all_affected_peaks.end());
      std::string peaks_str;
      size_t j = 0;
      for (const auto& peak : unique_peaks) {
        if (j > 0) peaks_str += ",";
        peaks_str += peak;
        j++;
      }
      associated_peaks_vec[i] = peaks_str;
    }

    // Set cascade_peak_genes and collect link mechanisms for cascade links
    if (!cascade_peak_gene_pairs.empty()) {
      // Remove duplicates
      std::unordered_set<std::string> unique_pairs(
          cascade_peak_gene_pairs.begin(), cascade_peak_gene_pairs.end());
      std::string pairs_str;
      size_t j = 0;
      for (const auto& pair : unique_pairs) {
        if (j > 0) pairs_str += ",";
        pairs_str += pair;
        j++;
      }
      cascade_peak_genes_vec[i] = pairs_str;

      // Collect mechanisms for cascade peak-gene links
      if (has_mechanism_col) {
        std::set<std::string> unique_mechanisms;
        for (const auto& pair : unique_pairs) {
          // pair is "peak->gene", convert to "peak::gene" for lookup
          size_t arrow_pos = pair.find("->");
          if (arrow_pos != std::string::npos) {
            std::string key =
                pair.substr(0, arrow_pos) + "::" + pair.substr(arrow_pos + 2);
            auto it = peak_gene_to_mechanism.find(key);
            if (it != peak_gene_to_mechanism.end()) {
              unique_mechanisms.insert(it->second);
            }
          }
        }
        if (!unique_mechanisms.empty()) {
          std::string mech_str;
          size_t k = 0;
          for (const auto& mech : unique_mechanisms) {
            if (k > 0) mech_str += ",";
            mech_str += mech;
            k++;
          }
          link_mechanism_vec[i] = mech_str;
        }
      }
    }

    // Map pattern to mechanism category using passed mapping
    if (pattern >= 1 && pattern <= pattern_to_mechanism.size()) {
      int mechanism_idx =
          pattern_to_mechanism[pattern - 1] - 1;  // Convert to 0-based index
      if (mechanism_idx >= 0 &&
          mechanism_idx < qtl_mechanism_categories.size()) {
        qtl_categories[i] = qtl_mechanism_categories[mechanism_idx];
      } else {
        qtl_categories[i] =
            qtl_mechanism_categories[7];  // Fallback: "No molQTL" (0-based
                                          // index 7 = 8th entry)
      }
    } else {
      qtl_categories[i] =
          qtl_mechanism_categories[7];  // Default: "No Detected QTL"
    }
  }

  return List::create(Named("variant_id") = variant_ids,
                      Named("qtl_pattern_number") = qtl_patterns,
                      Named("qtl_mechanism_category") = qtl_categories,
                      Named("associated_genes") = associated_genes_vec,
                      Named("associated_peaks") = associated_peaks_vec,
                      Named("cascade_peak_genes") = cascade_peak_genes_vec,
                      Named("peak_overlap") = peak_overlap,
                      Named("caqtl") = caqtl_status,
                      Named("peak_gene_link") = peak_gene_link_status,
                      Named("eqtl") = eqtl_status,
                      Named("link_mechanism") = link_mechanism_vec);
}

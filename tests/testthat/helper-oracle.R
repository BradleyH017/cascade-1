# helper-oracle.R
#
# Standalone reference oracle for CASCADE categorization logic.
# This file implements the categorization rules as plain R, without importing
# any cascade package functions. It serves as a ground-truth reference for
# testing that the package's C++ and R implementations produce correct results.
#
# NO library() calls — base R only.

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ORACLE_PBMC <- "predicted.celltype.l1.PBMC"

ORACLE_MYELOID <- c("predicted.celltype.l1.Mono", "predicted.celltype.l1.DC")

ORACLE_LYMPHOID <- c(
  "predicted.celltype.l1.NK",
  "predicted.celltype.l1.B",
  "predicted.celltype.l1.CD4_T",
  "predicted.celltype.l1.CD8_T",
  "predicted.celltype.l1.other_T"
)

ORACLE_TCELL <- c(
  "predicted.celltype.l1.CD4_T",
  "predicted.celltype.l1.CD8_T",
  "predicted.celltype.l1.other_T"
)

ORACLE_CATEGORIES <- c(
  "Cross-lineage shared",
  "Likely shared but underpowered",
  "Lineage-specific",
  "T-cell-specific",
  "Single cell-type",
  "No significance"
)

ORACLE_L2_TO_L1 <- list(
  "predicted.celltype.l2.CD14_Mono" = "predicted.celltype.l1.Mono",
  "predicted.celltype.l2.CD16_Mono" = "predicted.celltype.l1.Mono",
  "predicted.celltype.l2.cDC2" = "predicted.celltype.l1.DC",
  "predicted.celltype.l2.pDC" = "predicted.celltype.l1.DC",
  "predicted.celltype.l2.B_naive" = "predicted.celltype.l1.B",
  "predicted.celltype.l2.B_memory" = "predicted.celltype.l1.B",
  "predicted.celltype.l2.B_intermediate" = "predicted.celltype.l1.B",
  "predicted.celltype.l2.Plasmablast" = "predicted.celltype.l1.B",
  "predicted.celltype.l2.NK" = "predicted.celltype.l1.NK",
  "predicted.celltype.l2.NK_CD56bright" = "predicted.celltype.l1.NK",
  "predicted.celltype.l2.NK_Proliferating" = "predicted.celltype.l1.NK",
  "predicted.celltype.l2.CD4_Naive" = "predicted.celltype.l1.CD4_T",
  "predicted.celltype.l2.CD4_TCM" = "predicted.celltype.l1.CD4_T",
  "predicted.celltype.l2.CD4_TEM" = "predicted.celltype.l1.CD4_T",
  "predicted.celltype.l2.CD4_CTL" = "predicted.celltype.l1.CD4_T",
  "predicted.celltype.l2.Treg" = "predicted.celltype.l1.CD4_T",
  "predicted.celltype.l2.CD8_Naive" = "predicted.celltype.l1.CD8_T",
  "predicted.celltype.l2.CD8_TEM" = "predicted.celltype.l1.CD8_T",
  "predicted.celltype.l2.MAIT" = "predicted.celltype.l1.other_T",
  "predicted.celltype.l2.dnT" = "predicted.celltype.l1.other_T",
  "predicted.celltype.l2.gdT" = "predicted.celltype.l1.other_T",
  "predicted.celltype.l2.HSPC" = "predicted.celltype.l1.other",
  "predicted.celltype.l2.ILC" = "predicted.celltype.l1.other",
  "predicted.celltype.l2.Platelet" = "predicted.celltype.l1.other"
)

# ---------------------------------------------------------------------------
# oracle_map_to_l1
#
# Maps L2 cell type names to L1. If a CT is not found in the mapping, it is
# returned as-is (assumed to already be L1). Returns unique L1 CTs.
# ---------------------------------------------------------------------------

oracle_map_to_l1 <- function(cell_types, l2_to_l1 = ORACLE_L2_TO_L1) {
  mapped <- vapply(cell_types, function(ct) {
    if (ct %in% names(l2_to_l1)) {
      l2_to_l1[[ct]]
    } else {
      ct
    }
  }, character(1), USE.NAMES = FALSE)
  unique(mapped)
}

# ---------------------------------------------------------------------------
# oracle_categorize_specificity
#
# Pure-R reimplementation of categorize_single_feature_internal from
# src/categorization_core.cpp. Arguments:
#   sig_l1_cts    - character vector of significant L1 cell types
#   lfsr_values   - named numeric vector (names = CT names, values = LFSR);
#                   used for gray-zone checks
#   tested_l1_cts - if non-NULL, restricts which CTs are eligible for
#                   gray-zone checking (intersected with lfsr_values names)
#   lfsr_sig      - significance threshold (default 0.05)
#   lfsr_null     - null threshold for gray zone upper bound (default 0.5)
#
# Returns a single character string: one of ORACLE_CATEGORIES.
# ---------------------------------------------------------------------------

oracle_categorize_specificity <- function(sig_l1_cts,
                                          lfsr_values = NULL,
                                          tested_l1_cts = NULL,
                                          lfsr_sig = 0.05,
                                          lfsr_null = 0.5) {
  # Step 1: separate PBMC from other significant CTs
  # Use filtering (not setdiff) to preserve multiplicity, matching C++ behavior
  sig_non_pbmc <- sig_l1_cts[sig_l1_cts != ORACLE_PBMC]
  has_pbmc <- ORACLE_PBMC %in% sig_l1_cts

  # Step 2: only PBMC significant
  if (length(sig_non_pbmc) == 0 && has_pbmc) {
    return("Likely shared but underpowered")
  }

  # Step 3: no significant CTs at all
  if (length(sig_non_pbmc) == 0) {
    return("No significance")
  }

  # Step 4: count lineage membership among significant non-PBMC CTs
  myeloid_sig <- intersect(sig_non_pbmc, ORACLE_MYELOID)
  lymphoid_sig <- intersect(sig_non_pbmc, ORACLE_LYMPHOID)
  tcell_sig <- intersect(sig_non_pbmc, ORACLE_TCELL)

  has_myeloid <- length(myeloid_sig) > 0
  has_lymphoid <- length(lymphoid_sig) > 0

  # Step 5: both lineages present
  if (has_myeloid && has_lymphoid) {
    return("Cross-lineage shared")
  }

  # Helper: check if any CT in check_cts has gray-zone LFSR
  has_gray_zone <- function(check_cts) {
    if (is.null(lfsr_values) || length(lfsr_values) == 0) {
      return(FALSE)
    }
    # Optionally restrict to tested CTs
    if (!is.null(tested_l1_cts)) {
      check_cts <- intersect(check_cts, tested_l1_cts)
    }
    for (ct in check_cts) {
      if (ct %in% names(lfsr_values)) {
        val <- lfsr_values[[ct]]
        # Handle NaN/NA: in C++, NaN comparisons return false, so NaN is not gray zone
        if (!is.na(val) && !is.nan(val) && val >= lfsr_sig && val < lfsr_null) {
          return(TRUE)
        }
      }
    }
    FALSE
  }

  # Step 6: exactly 1 significant CT (excluding PBMC)
  if (length(sig_non_pbmc) == 1) {
    # Check ALL non-sig, non-PBMC CTs that appear in lfsr_values
    # "other" CTs (not in myeloid/lymphoid) ARE checked here.
    all_lfsr_cts <- names(lfsr_values)
    nonsig_nonpbmc <- setdiff(all_lfsr_cts, c(sig_non_pbmc, ORACLE_PBMC))
    if (!is.null(tested_l1_cts)) {
      nonsig_nonpbmc <- intersect(nonsig_nonpbmc, tested_l1_cts)
    }
    if (has_gray_zone(nonsig_nonpbmc)) {
      return("Likely shared but underpowered")
    }
    return("Single cell-type")
  }

  # Step 7: T-cells only (all lymphoid_sig are T-cells, no myeloid)
  if (!has_myeloid && has_lymphoid &&
    length(lymphoid_sig) == length(tcell_sig) &&
    all(lymphoid_sig %in% ORACLE_TCELL)) {
    # Check CTs in myeloid OR (lymphoid AND NOT T-cell) — "is_other_lineage"
    # Note: "other" CTs (not in myeloid/lymphoid) are NOT checked here.
    other_lineage_cts <- c(ORACLE_MYELOID, setdiff(ORACLE_LYMPHOID, ORACLE_TCELL))
    if (has_gray_zone(other_lineage_cts)) {
      return("Likely shared but underpowered")
    }
    return("T-cell-specific")
  }

  # Step 8: one lineage only (myeloid XOR lymphoid)
  if (has_myeloid && !has_lymphoid) {
    # Check LFSR for CTs in the ABSENT lineage (lymphoid) only
    # "other" CTs are NOT checked.
    if (has_gray_zone(ORACLE_LYMPHOID)) {
      return("Likely shared but underpowered")
    }
    return("Lineage-specific")
  }
  if (has_lymphoid && !has_myeloid) {
    # Check LFSR for CTs in the ABSENT lineage (myeloid) only
    # "other" CTs are NOT checked.
    if (has_gray_zone(ORACLE_MYELOID)) {
      return("Likely shared but underpowered")
    }
    return("Lineage-specific")
  }

  # Step 9: default fallback
  "No significance"
}

# ---------------------------------------------------------------------------
# oracle_derive_pattern
#
# Implements the 25-pattern QTL mechanism decision tree.
# All arguments are logical scalars.
#
# Returns an integer pattern number (1-25).
# ---------------------------------------------------------------------------

oracle_derive_pattern <- function(has_overlap,
                                  has_caqtl_overlap,
                                  has_caqtl_nonoverlap,
                                  has_link_overlap,
                                  has_link_nonoverlap,
                                  has_eqtl,
                                  eqtl_for_linked_overlap,
                                  eqtl_for_linked_nonoverlap) {
  if (has_overlap) {
    if (has_caqtl_overlap) {
      if (has_link_overlap) {
        if (has_eqtl && eqtl_for_linked_overlap) {
          return(1L)
        } # Local Cascade
        if (has_eqtl && !eqtl_for_linked_overlap) {
          return(6L)
        } # caQTL+eQTL (overlap link, discordant)
        if (!has_eqtl) {
          return(13L)
        } # Only caQTL (overlap link)
      } else {
        if (has_eqtl) {
          return(10L)
        } # caQTL+eQTL (no link)
        if (!has_eqtl) {
          return(17L)
        } # Only caQTL (overlap, no link)
      }
    } else if (has_caqtl_nonoverlap) {
      if (has_link_overlap) {
        if (has_eqtl && eqtl_for_linked_overlap) {
          return(2L)
        } # Positional Cascade (overlap link, non-overlap caQTL)
        if (has_eqtl && !eqtl_for_linked_overlap) {
          return(7L)
        } # caQTL+eQTL (overlap link, non-overlap caQTL, discordant)
        if (!has_eqtl) {
          return(14L)
        } # Only caQTL (overlap link, non-overlap caQTL)
      } else if (has_link_nonoverlap) {
        if (has_eqtl && eqtl_for_linked_nonoverlap) {
          return(4L)
        } # Positional Cascade (non-overlap link)
        if (has_eqtl && !eqtl_for_linked_nonoverlap) {
          return(8L)
        } # caQTL+eQTL (non-overlap link, discordant)
        if (!has_eqtl) {
          return(15L)
        } # Only caQTL (non-overlap link)
      } else {
        if (has_eqtl) {
          return(11L)
        } # caQTL+eQTL (non-overlap caQTL, no link)
        if (!has_eqtl) {
          return(18L)
        } # Only caQTL (non-overlap caQTL, no link)
      }
    } else {
      # no caQTL, but overlap
      if (has_link_overlap) {
        if (has_eqtl && eqtl_for_linked_overlap) {
          return(3L)
        } # Positional Cascade (overlap link)
        if (has_eqtl && !eqtl_for_linked_overlap) {
          return(20L)
        } # Only eQTL (overlap link)
        if (!has_eqtl) {
          return(23L)
        } # No molQTL (overlap link)
      } else {
        if (has_eqtl) {
          return(21L)
        } # Only eQTL (no link)
        if (!has_eqtl) {
          return(24L)
        } # No molQTL (no link)
      }
    }
  } else {
    # no overlap
    if (has_caqtl_nonoverlap) {
      if (has_link_nonoverlap) {
        if (has_eqtl && eqtl_for_linked_nonoverlap) {
          return(5L)
        } # Distal Cascade
        if (has_eqtl && !eqtl_for_linked_nonoverlap) {
          return(9L)
        } # caQTL+eQTL (no overlap, discordant)
        if (!has_eqtl) {
          return(16L)
        } # Only caQTL (no overlap, link)
      } else {
        if (has_eqtl) {
          return(12L)
        } # caQTL+eQTL (no overlap, no link)
        if (!has_eqtl) {
          return(19L)
        } # Only caQTL (no overlap, no link)
      }
    } else {
      # no caQTL, no overlap
      if (has_eqtl) {
        return(22L)
      } # Only eQTL (no overlap)
      if (!has_eqtl) {
        return(25L)
      } # No molQTL (no overlap)
    }
  }
}

# ---------------------------------------------------------------------------
# oracle_pattern_to_mechanism
#
# Maps a pattern number (1-25) to its mechanism string.
# ---------------------------------------------------------------------------

oracle_pattern_to_mechanism <- function(pattern_number) {
  if (pattern_number == 1L) {
    return("Local Cascade")
  }
  if (pattern_number >= 2L && pattern_number <= 4L) {
    return("Positional Cascade")
  }
  if (pattern_number == 5L) {
    return("Distal Cascade")
  }
  if (pattern_number >= 6L && pattern_number <= 12L) {
    return("caQTL + eQTL (No Link)")
  }
  if (pattern_number >= 13L && pattern_number <= 16L) {
    return("Only caQTL (With Link)")
  }
  if (pattern_number >= 17L && pattern_number <= 19L) {
    return("Only caQTL (No Link)")
  }
  if (pattern_number >= 20L && pattern_number <= 22L) {
    return("Only eQTL")
  }
  if (pattern_number >= 23L && pattern_number <= 25L) {
    return("No molQTL")
  }
  stop("Invalid pattern number: ", pattern_number)
}

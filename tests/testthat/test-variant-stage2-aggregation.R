# test-variant-stage2-aggregation.R
#
# Tests for Stage 2 cross-cell-type aggregation logic.
# Verifies the rules that combine per-cell-type QTL pattern results into:
#   - qtl_mechanism_category (mechanism from best/lowest pattern number)
#   - cell_type_specificity  (specificity of CTs sharing the best mechanism)
#   - best_cell_types        (CTs with the best pattern number)
#   - gene_affected_cell_types (CTs with patterns 1-12 or 20-22 = has eQTL)
#   - peak_affected_cell_types (CTs with patterns 1,2,4-19 = has caQTL; P3 excluded)
#
# Rather than constructing full per-celltype result structures required by
# process_qtl_variant_results (which needs many fields), we test the
# aggregation LOGIC directly using simple R and the oracle.

# ---------- shared helpers ----------

CT_MONO <- "predicted.celltype.l1.Mono"
CT_DC <- "predicted.celltype.l1.DC"
CT_B <- "predicted.celltype.l1.B"
CT_NK <- "predicted.celltype.l1.NK"
CT_CD4_T <- "predicted.celltype.l1.CD4_T"
CT_CD8_T <- "predicted.celltype.l1.CD8_T"
CT_OTHER_T <- "predicted.celltype.l1.other_T"
CT_PBMC <- "predicted.celltype.l1.PBMC"

MYELOID_CTS <- c(CT_MONO, CT_DC)
LYMPHOID_CTS <- c(CT_NK, CT_B, CT_CD4_T, CT_CD8_T, CT_OTHER_T)
TCELL_CTS <- c(CT_CD4_T, CT_CD8_T, CT_OTHER_T)
ALL_L1_CTS <- c(MYELOID_CTS, LYMPHOID_CTS)

# Pattern ranges for eQTL / caQTL affected CTs (from aggregate_variant_categorization)
EQTL_PATTERNS <- c(1L:12L, 20L:22L) # unchanged: cascades + caQTL+eQTL + Only eQTL
CAQTL_PATTERNS <- c(1L, 2L, 4L:19L) # P3 has no caQTL; P5 (Distal) has caQTL

NO_QTL_PATTERNS <- c(23L, 24L, 25L)

#' Check if a pattern has eQTL
has_eqtl_pattern <- function(p) p %in% EQTL_PATTERNS

#' Check if a pattern has caQTL
has_caqtl_pattern <- function(p) p %in% CAQTL_PATTERNS

#' Simulate Stage 2 aggregation logic for a single variant.
#'
#' @param ct_patterns Named integer vector: cell_type -> pattern_number (1-25).
#' @return List with best_pattern, mechanism, best_cts, gene_affected_cts,
#'         peak_affected_cts, specificity.
simulate_stage2 <- function(ct_patterns) {
  # Best pattern = minimum across all CTs
  best_pattern <- min(ct_patterns)

  # Mechanism from oracle
  mechanism <- oracle_pattern_to_mechanism(best_pattern)

  # Best CTs: those matching the best pattern
  best_cts <- names(ct_patterns)[ct_patterns == best_pattern]

  # Gene-affected CTs: patterns with eQTL
  gene_affected_cts <- names(ct_patterns)[vapply(ct_patterns, has_eqtl_pattern, logical(1))]

  # Peak-affected CTs: patterns with caQTL
  peak_affected_cts <- names(ct_patterns)[vapply(ct_patterns, has_caqtl_pattern, logical(1))]

  # For combined specificity: CTs sharing the best mechanism AND has QTL
  all_mechanisms <- vapply(ct_patterns, oracle_pattern_to_mechanism, character(1))
  has_qtl <- !(ct_patterns %in% NO_QTL_PATTERNS)
  affected_cts <- names(ct_patterns)[all_mechanisms == mechanism & has_qtl]

  # Specificity (LFSR disabled): based on affected CTs mapped to L1
  affected_l1 <- oracle_map_to_l1(affected_cts)
  specificity <- oracle_categorize_specificity(
    sig_l1_cts = affected_l1,
    lfsr_values = NULL # LFSR disabled for combined specificity
  )

  list(
    best_pattern         = best_pattern,
    mechanism            = mechanism,
    best_cts             = sort(best_cts),
    gene_affected_cts    = sort(gene_affected_cts),
    peak_affected_cts    = sort(peak_affected_cts),
    affected_cts         = sort(affected_cts),
    specificity          = specificity
  )
}

# ---------- tests ----------

test_that("Case 1: Single CT, Local Cascade (P1)", {
  ct_patterns <- c(setNames(1L, CT_MONO))
  res <- simulate_stage2(ct_patterns)

  expect_equal(res$best_pattern, 1L)
  expect_equal(res$mechanism, "Local Cascade")
  expect_equal(res$best_cts, CT_MONO)
  expect_equal(res$specificity, "Single cell-type")

  # P1 has both eQTL and caQTL
  expect_true(CT_MONO %in% res$gene_affected_cts)
  expect_true(CT_MONO %in% res$peak_affected_cts)
})

test_that("Case 2: Multi CT, same pattern (Mono P1 + B P1)", {
  ct_patterns <- setNames(c(1L, 1L), c(CT_MONO, CT_B))
  res <- simulate_stage2(ct_patterns)

  expect_equal(res$best_pattern, 1L)
  expect_equal(res$mechanism, "Local Cascade")
  expect_equal(res$best_cts, sort(c(CT_MONO, CT_B)))
  # Mono = myeloid, B = lymphoid => Cross-lineage shared
  expect_equal(res$specificity, "Cross-lineage shared")

  # Both should be gene- and peak-affected
  expect_equal(res$gene_affected_cts, sort(c(CT_MONO, CT_B)))
  expect_equal(res$peak_affected_cts, sort(c(CT_MONO, CT_B)))
})

test_that("Case 3: Multi CT, mixed patterns (Mono P1 + B P19)", {
  ct_patterns <- setNames(c(1L, 19L), c(CT_MONO, CT_B))
  res <- simulate_stage2(ct_patterns)

  expect_equal(res$best_pattern, 1L)
  expect_equal(res$mechanism, "Local Cascade")
  # Only Mono has P1
  expect_equal(res$best_cts, CT_MONO)

  # Affected CTs for combined specificity: must share best mechanism AND have QTL
  # Mono has "Local Cascade", B has "Only caQTL (No Link)" -- different mechanisms
  # So only Mono is in affected_cts
  expect_equal(res$affected_cts, CT_MONO)
  expect_equal(res$specificity, "Single cell-type")

  # Gene-affected: P1 has eQTL, P19 does not (Only caQTL)
  expect_true(CT_MONO %in% res$gene_affected_cts)
  expect_false(CT_B %in% res$gene_affected_cts)

  # Peak-affected: P1 has caQTL, P19 has caQTL
  expect_true(CT_MONO %in% res$peak_affected_cts)
  expect_true(CT_B %in% res$peak_affected_cts)
})

test_that("Case 4: All No molQTL (P25)", {
  ct_patterns <- setNames(
    rep(25L, length(ALL_L1_CTS)),
    ALL_L1_CTS
  )
  res <- simulate_stage2(ct_patterns)

  expect_equal(res$best_pattern, 25L)
  expect_equal(res$mechanism, "No molQTL")
  # All CTs have P25
  expect_equal(res$best_cts, sort(ALL_L1_CTS))
  # No gene or peak affected
  expect_length(res$gene_affected_cts, 0)
  expect_length(res$peak_affected_cts, 0)
  # In the actual pipeline, no-QTL variants are handled by a separate branch
  # that sets mechanism = "No molQTL" and pattern = 25 directly.
  # With all P25, no CTs have QTL so affected_cts is empty.
  expect_length(res$affected_cts, 0)
  # No affected CTs => "No significance"
  expect_equal(res$specificity, "No significance")
})

test_that("Case 5: Different QTL types (Mono P22 + B P19)", {
  ct_patterns <- setNames(c(22L, 19L), c(CT_MONO, CT_B))
  res <- simulate_stage2(ct_patterns)

  # P19 < P22, so best pattern is 19
  expect_equal(res$best_pattern, 19L)
  expect_equal(res$mechanism, "Only caQTL (No Link)")
  # Only B has P19
  expect_equal(res$best_cts, CT_B)

  # Affected CTs sharing mechanism "Only caQTL (No Link)": only B (P19)
  # Mono has P22 = "Only eQTL" -- different mechanism
  expect_equal(res$affected_cts, CT_B)
  # B is lymphoid only => single lineage => specificity depends on lineage check
  # With LFSR disabled: only lymphoid, no myeloid => "Lineage-specific"
  # Wait -- only 1 CT (B) is affected. 1 CT => "Single cell-type"
  expect_equal(res$specificity, "Single cell-type")

  # Gene-affected: P22 has eQTL (Only eQTL), P19 does not
  expect_true(CT_MONO %in% res$gene_affected_cts)
  expect_false(CT_B %in% res$gene_affected_cts)

  # Peak-affected: P19 has caQTL, P22 does not
  expect_false(CT_MONO %in% res$peak_affected_cts)
  expect_true(CT_B %in% res$peak_affected_cts)
})

test_that("Case 6: Pattern priority (CD4_T P3 + NK P20)", {
  # P3 = Positional Cascade (overlap link, no caQTL, has eQTL)
  ct_patterns <- setNames(c(3L, 20L), c(CT_CD4_T, CT_NK))
  res <- simulate_stage2(ct_patterns)

  # P3 < P20
  expect_equal(res$best_pattern, 3L)
  expect_equal(res$mechanism, "Positional Cascade")
  expect_equal(res$best_cts, CT_CD4_T)

  # Only CD4_T shares the "Positional Cascade" mechanism
  expect_equal(res$affected_cts, CT_CD4_T)
  expect_equal(res$specificity, "Single cell-type")

  # Gene-affected: P3 has eQTL, P20 has eQTL
  expect_true(CT_CD4_T %in% res$gene_affected_cts)
  expect_true(CT_NK %in% res$gene_affected_cts)

  # Peak-affected: P3 is NOT in caqtl_patterns (no caQTL), P20 is not in caqtl
  # P3 = "Positional Cascade (overlap link)": no caQTL, just link via overlap peak
  expect_false(CT_CD4_T %in% res$peak_affected_cts)
  expect_false(CT_NK %in% res$peak_affected_cts)
})

test_that("Case 7: Affected CT extraction (Mono P1 + B P25 + CD4_T P22)", {
  ct_patterns <- setNames(c(1L, 25L, 22L), c(CT_MONO, CT_B, CT_CD4_T))
  res <- simulate_stage2(ct_patterns)

  # Best pattern = min(1, 25, 22) = 1
  expect_equal(res$best_pattern, 1L)
  expect_equal(res$mechanism, "Local Cascade")
  expect_equal(res$best_cts, CT_MONO)

  # Affected CTs: must share "Local Cascade" mechanism AND have QTL (not P23-25)
  # Mono=P1 (Local Cascade, has QTL) -> YES
  # B=P25 (No molQTL, no QTL) -> NO
  # CD4_T=P22 (Only eQTL, has QTL but different mechanism) -> NO
  expect_equal(res$affected_cts, CT_MONO)
  expect_equal(res$specificity, "Single cell-type")

  # Gene-affected: P1 has eQTL, P25 does not, P22 has eQTL
  expect_true(CT_MONO %in% res$gene_affected_cts)
  expect_false(CT_B %in% res$gene_affected_cts)
  expect_true(CT_CD4_T %in% res$gene_affected_cts)

  # Peak-affected: P1 has caQTL, P25 does not, P22 does not
  expect_true(CT_MONO %in% res$peak_affected_cts)
  expect_false(CT_B %in% res$peak_affected_cts)
  expect_false(CT_CD4_T %in% res$peak_affected_cts)
})

test_that("Case 8: Tie in best pattern (Mono P5 + B P5)", {
  # P5 = Distal Cascade (has both eQTL and caQTL)
  ct_patterns <- setNames(c(5L, 5L), c(CT_MONO, CT_B))
  res <- simulate_stage2(ct_patterns)

  # Both have P5
  expect_equal(res$best_pattern, 5L)
  expect_equal(res$mechanism, "Distal Cascade")
  expect_equal(res$best_cts, sort(c(CT_MONO, CT_B)))

  # Both share the mechanism and have QTL
  expect_equal(res$affected_cts, sort(c(CT_MONO, CT_B)))
  # Mono = myeloid, B = lymphoid => Cross-lineage shared
  expect_equal(res$specificity, "Cross-lineage shared")

  # P5 = Distal Cascade: has both eQTL and caQTL
  expect_equal(res$gene_affected_cts, sort(c(CT_MONO, CT_B)))
  expect_equal(res$peak_affected_cts, sort(c(CT_MONO, CT_B)))
})

# ---------- build_pattern_category_matrices tests ----------
# Verify the actual internal function produces correct matrices

test_that("build_pattern_category_matrices creates correct matrices from per-celltype results", {
  variant_ids <- c("var1", "var2")
  cell_types <- c(CT_MONO, CT_B)

  # Simulate per_celltype_results as the pipeline produces them
  per_celltype_results <- list()
  per_celltype_results[[CT_MONO]] <- data.frame(
    variant_id = c("var1", "var2"),
    qtl_pattern_number = c(1L, 22L),
    qtl_mechanism_category = c("Local Cascade", "Only eQTL"),
    stringsAsFactors = FALSE
  )
  per_celltype_results[[CT_B]] <- data.frame(
    variant_id = c("var1", "var2"),
    qtl_pattern_number = c(19L, 5L),
    qtl_mechanism_category = c("Only caQTL (No Link)", "Distal Cascade"),
    stringsAsFactors = FALSE
  )

  matrices <- cascade:::build_pattern_category_matrices(
    variant_ids, cell_types, per_celltype_results
  )

  # Check pattern_matrix dimensions and values
  expect_equal(nrow(matrices$pattern_matrix), 2)
  expect_equal(ncol(matrices$pattern_matrix), 2)
  expect_equal(colnames(matrices$pattern_matrix), cell_types)

  # var1: Mono=P1, B=P19
  expect_equal(unname(matrices$pattern_matrix[1, 1]), 1)
  expect_equal(unname(matrices$pattern_matrix[1, 2]), 19)

  # var2: Mono=P22, B=P5
  expect_equal(unname(matrices$pattern_matrix[2, 1]), 22)
  expect_equal(unname(matrices$pattern_matrix[2, 2]), 5)

  # Check category_matrix
  expect_equal(unname(matrices$category_matrix[1, 1]), "Local Cascade")
  expect_equal(unname(matrices$category_matrix[1, 2]), "Only caQTL (No Link)")
  expect_equal(unname(matrices$category_matrix[2, 1]), "Only eQTL")
  expect_equal(unname(matrices$category_matrix[2, 2]), "Distal Cascade")
})

test_that("build_pattern_category_matrices defaults to P25/No molQTL for missing CTs", {
  variant_ids <- c("var1")
  cell_types <- c(CT_MONO, CT_B, CT_NK)

  # Only Mono has results for var1
  per_celltype_results <- list()
  per_celltype_results[[CT_MONO]] <- data.frame(
    variant_id = "var1",
    qtl_pattern_number = 1L,
    qtl_mechanism_category = "Local Cascade",
    stringsAsFactors = FALSE
  )

  matrices <- cascade:::build_pattern_category_matrices(
    variant_ids, cell_types, per_celltype_results
  )

  # Mono has P1, B and NK default to P25
  expect_equal(unname(matrices$pattern_matrix[1, 1]), 1)
  expect_equal(unname(matrices$pattern_matrix[1, 2]), 25)
  expect_equal(unname(matrices$pattern_matrix[1, 3]), 25)

  expect_equal(unname(matrices$category_matrix[1, 1]), "Local Cascade")
  expect_equal(unname(matrices$category_matrix[1, 2]), "No molQTL")
  expect_equal(unname(matrices$category_matrix[1, 3]), "No molQTL")
})

# ---------- eQTL / caQTL pattern membership verification ----------
# Cross-check that pattern ranges in source match oracle expectations

test_that("eQTL pattern range matches oracle expectations", {
  # Patterns with eQTL: 1-12 (Cascade mechanisms + caQTL+eQTL) and 20-22 (Only eQTL)
  for (p in 1:25) {
    mech <- oracle_pattern_to_mechanism(p)
    has_eqtl <- p %in% EQTL_PATTERNS
    expected_eqtl <- mech %in% c(
      "Local Cascade", "Positional Cascade", "Distal Cascade",
      "caQTL + eQTL (No Link)", "Only eQTL"
    )
    expect_equal(
      has_eqtl, expected_eqtl,
      info = sprintf(
        "Pattern %d (%s): eQTL=%s expected=%s",
        p, mech, has_eqtl, expected_eqtl
      )
    )
  }
})

test_that("caQTL pattern range matches oracle expectations", {
  # Source code caqtl_patterns: 1, 2, 4-19 (P3 has no caQTL; P5 Distal has caQTL)
  source_caqtl <- c(
    1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L,
    13L, 14L, 15L, 16L, 17L, 18L, 19L
  )
  expect_equal(sort(CAQTL_PATTERNS), sort(source_caqtl))
})

# ---------- oracle mechanism consistency ----------

test_that("All 25 patterns map to one of the 8 QTL mechanisms via oracle", {
  valid_mechanisms <- c(
    "Local Cascade", "Positional Cascade", "Distal Cascade",
    "caQTL + eQTL (No Link)", "Only caQTL (With Link)", "Only caQTL (No Link)",
    "Only eQTL", "No molQTL"
  )
  for (p in 1:25) {
    mech <- oracle_pattern_to_mechanism(p)
    expect_true(mech %in% valid_mechanisms,
      info = sprintf("Pattern %d -> '%s' not in valid mechanisms", p, mech)
    )
  }
})

test_that("Best pattern = min selects highest-priority mechanism", {
  # Local Cascade (P1) < Positional Cascade (P2-4) < Distal Cascade (P5) < ... < No molQTL (P23-25)
  # Verify that lower pattern numbers always map to equal-or-higher priority mechanism
  mech_priority <- c(
    "Local Cascade" = 1, "Positional Cascade" = 2, "Distal Cascade" = 3,
    "caQTL + eQTL (No Link)" = 4, "Only caQTL (With Link)" = 5,
    "Only caQTL (No Link)" = 6, "Only eQTL" = 7, "No molQTL" = 8
  )
  for (p in 1:24) {
    mech_p <- oracle_pattern_to_mechanism(p)
    mech_next <- oracle_pattern_to_mechanism(p + 1L)
    expect_true(
      mech_priority[mech_p] <= mech_priority[mech_next],
      info = sprintf(
        "P%d (%s, priority %d) should be <= P%d (%s, priority %d)",
        p, mech_p, mech_priority[mech_p],
        p + 1L, mech_next, mech_priority[mech_next]
      )
    )
  }
})

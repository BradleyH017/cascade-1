# test-feature-categorization-logic.R
#
# Tests for categorize_features_from_acat() C++ function.
# Each test constructs a minimal acat_matrix data.frame and verifies
# the ACAT q-value -> specificity category mapping.

# ---------- shared helpers ----------

# L1 cell type column names
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

ALL_L1_CTS <- c(MYELOID_CTS, LYMPHOID_CTS) # excludes PBMC (handled separately)

#' Build a minimal acat_matrix data.frame.
#' @param feature_id Character scalar.
#' @param sig_values Named list of cell type -> q-value for significant entries.
#' @param default_q  Default q-value for non-significant cell types.
#' @param ct_cols    Character vector of cell type column names to include.
#' @param include_pbmc Logical; if TRUE add a PBMC column.
#' @param pbmc_q     PBMC q-value.
make_acat <- function(feature_id,
                      sig_values = list(),
                      default_q = 0.5,
                      ct_cols = ALL_L1_CTS,
                      include_pbmc = FALSE,
                      pbmc_q = 0.5) {
  df <- data.frame(feature_id = feature_id, stringsAsFactors = FALSE)
  for (ct in ct_cols) {
    df[[ct]] <- if (!is.null(sig_values[[ct]])) sig_values[[ct]] else default_q
  }
  if (include_pbmc) {
    df[[CT_PBMC]] <- pbmc_q
  }
  df
}

# ---------- tests ----------

test_that("Case 1 (G1): Cross-lineage shared - Myeloid + Lymphoid significant", {
  acat <- make_acat(
    "G1",
    sig_values = setNames(list(0.01, 0.02), c(CT_MONO, CT_B))
  )
  res <- categorize_features_from_acat(
    acat_matrix = acat,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC,
    other_cts = character(0),
    sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "Cross-lineage shared")
})

test_that("Case 2 (G2): Lineage-specific - Myeloid only", {
  acat <- make_acat(
    "G2",
    sig_values = setNames(list(0.01, 0.03), c(CT_MONO, CT_DC))
  )
  res <- categorize_features_from_acat(
    acat_matrix = acat,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC,
    other_cts = character(0),
    sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "Lineage-specific")
})

test_that("Case 3 (G3): No significance - all q > 0.05", {
  acat <- make_acat("G3")
  res <- categorize_features_from_acat(
    acat_matrix = acat,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC,
    other_cts = character(0),
    sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "No significance")
})

test_that("Case 4 (G4): Single cell-type - one T cell significant, no LFSR", {
  acat <- make_acat(
    "G4",
    sig_values = setNames(list(0.001), CT_CD4_T)
  )
  res <- categorize_features_from_acat(
    acat_matrix = acat,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC,
    other_cts = character(0),
    sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "Single cell-type")
})

test_that("Case 5 (G5): Likely shared but underpowered - PBMC only significant", {
  acat <- make_acat(
    "G5",
    include_pbmc = TRUE,
    pbmc_q       = 0.01
  )
  res <- categorize_features_from_acat(
    acat_matrix = acat,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC,
    other_cts = character(0),
    sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "Likely shared but underpowered")
})

test_that("Case 6 (G6): Single cell-type - L2 subtypes map to same L1 (Mono)", {
  # Use L2 column names for Mono subtypes; all other L1 CTs are non-significant.
  l2_cd14 <- "predicted.celltype.l2.CD14_Mono"
  l2_cd16 <- "predicted.celltype.l2.CD16_Mono"

  ct_cols_l2 <- c(l2_cd14, l2_cd16, CT_B, CT_NK, CT_CD4_T, CT_CD8_T, CT_OTHER_T)

  acat <- make_acat(
    "G6",
    sig_values = setNames(list(0.01, 0.02), c(l2_cd14, l2_cd16)),
    ct_cols    = ct_cols_l2
  )

  l2_map <- list(
    "predicted.celltype.l2.CD14_Mono" = "predicted.celltype.l1.Mono",
    "predicted.celltype.l2.CD16_Mono" = "predicted.celltype.l1.Mono"
  )

  res <- categorize_features_from_acat(
    acat_matrix = acat,
    l2_to_l1_mapping = l2_map,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC,
    other_cts = character(0),
    sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "Single cell-type")
})

test_that("Case 7 (G7): Boundary - q=0.049 is significant -> Single cell-type", {
  acat <- make_acat(
    "G7",
    sig_values = setNames(list(0.049), CT_MONO)
  )
  res <- categorize_features_from_acat(
    acat_matrix = acat,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC,
    other_cts = character(0),
    sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "Single cell-type")
})

test_that("Case 8 (G8): Boundary - q=0.051 is not significant -> No significance", {
  acat <- make_acat(
    "G8",
    sig_values = setNames(list(0.051), CT_MONO)
  )
  res <- categorize_features_from_acat(
    acat_matrix = acat,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC,
    other_cts = character(0),
    sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "No significance")
})

# ---------- LFSR integration tests (new long-format interface) ----------

test_that("LFSR gray zone through C++ feature path: single CT + gray zone → underpowered", {
  acat <- make_acat("G_lfsr1", sig_values = setNames(list(0.01), CT_MONO))
  lfsr_lookup <- data.table::data.table(
    feature_id = "G_lfsr1",
    cell_type = CT_B,
    min_lfsr = 0.15
  )
  res <- categorize_features_from_acat(
    acat_matrix = acat, lfsr_lookup = lfsr_lookup,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC, other_cts = character(0), sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "Likely shared but underpowered")
})

test_that("LFSR null through C++ feature path: single CT + null LFSR → single cell-type", {
  acat <- make_acat("G_lfsr2", sig_values = setNames(list(0.01), CT_MONO))
  lfsr_lookup <- data.table::data.table(
    feature_id = "G_lfsr2",
    cell_type = CT_B,
    min_lfsr = 0.8
  )
  res <- categorize_features_from_acat(
    acat_matrix = acat, lfsr_lookup = lfsr_lookup,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC, other_cts = character(0), sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "Single cell-type")
})

test_that("LFSR at boundary through C++ feature path: LFSR=0.05 → gray zone", {
  acat <- make_acat("G_lfsr3", sig_values = setNames(list(0.01), CT_MONO))
  lfsr_lookup <- data.table::data.table(
    feature_id = "G_lfsr3",
    cell_type = CT_B,
    min_lfsr = 0.05
  )
  res <- categorize_features_from_acat(
    acat_matrix = acat, lfsr_lookup = lfsr_lookup,
    lineage_groups = list(MYELOID_CTS, LYMPHOID_CTS),
    subgroup_levels = list(list(TCELL_CTS)),
    bulk_cts = CT_PBMC, other_cts = character(0), sig_threshold = 0.05,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(res$primary_category, "Likely shared but underpowered")
})

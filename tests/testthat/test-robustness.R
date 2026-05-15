# test-robustness.R
#
# Edge-case robustness tests for CASCADE categorization logic.

test_that("ACAT q-value all NA produces No significance", {
  # Build a minimal ACAT matrix where the only gene has NA q-values
  acat_matrix <- data.frame(
    feature_id = "gene1",
    predicted.celltype.l1.Mono = NA_real_,
    predicted.celltype.l1.NK = NA_real_,
    stringsAsFactors = FALSE
  )


  result <- categorize_features_from_acat(
    acat_matrix = acat_matrix,
    feature_type = "gene",
    lineage_groups = list(ORACLE_MYELOID, ORACLE_LYMPHOID),
    subgroup_levels = list(list(ORACLE_TCELL)),
    bulk_cts = ORACLE_PBMC,
    other_cts = character(0),
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )

  expect_equal(nrow(result), 1L)
  expect_equal(result$primary_category[1], "No significance")
  expect_equal(result$n_significant_cts[1], 0L)
})

test_that("LFSR NaN does not trigger gray zone", {
  # One significant CT (Mono), LFSR for NK is NaN.
  # In C++, NaN comparisons (>= and <) return false, so NaN should NOT
  # satisfy the gray zone condition. Result should be "Single cell-type".
  result <- categorize_cell_specificity(
    significant_cts = "predicted.celltype.l1.Mono",
    lfsr_values = NaN,
    lfsr_names = "predicted.celltype.l1.NK",
    lineage_groups = list(ORACLE_MYELOID, ORACLE_LYMPHOID),
    subgroup_levels = list(list(ORACLE_TCELL)),
    bulk_cts = ORACLE_PBMC,
    other_cts = character(0),
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )
  expect_equal(result, "Single cell-type")
})

test_that("empty ACAT data.table returns 0-row output with correct schema", {
  acat_matrix <- data.frame(
    feature_id = character(0),
    predicted.celltype.l1.Mono = numeric(0),
    predicted.celltype.l1.NK = numeric(0),
    stringsAsFactors = FALSE
  )

  result <- categorize_features_from_acat(
    acat_matrix = acat_matrix,
    feature_type = "gene",
    lineage_groups = list(ORACLE_MYELOID, ORACLE_LYMPHOID),
    subgroup_levels = list(list(ORACLE_TCELL)),
    bulk_cts = ORACLE_PBMC,
    other_cts = character(0),
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels
  )

  expect_equal(nrow(result), 0L)
  expected_cols <- c(
    "feature_id", "primary_category", "cell_type_specificity_pattern",
    "significant_cts", "tested_cts", "n_significant_cts", "n_tested_cts"
  )
  expect_equal(names(result), expected_cols)
})

test_that("unknown cell type counted as significant but belongs to no lineage", {
  # "predicted.celltype.l1.Alien" is not in any lineage set (myeloid/lymphoid/tcell).
  # It counts toward sig_non_pbmc (size=2), so we skip the single-CT branch.
  # Mono is myeloid, Alien is not — only myeloid lineage represented.
  # Result: "Lineage-specific" (myeloid only, lymphoid empty).
  result <- oracle_categorize_specificity(
    sig_l1_cts = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.Alien")
  )
  expect_equal(result, "Lineage-specific")
})

test_that("config hash changes when parameters change, not just thresholds", {
  config_a <- list(
    file_patterns = list(acat = "acat*.tsv"),
    cell_types = c("Mono", "NK"),
    parameters = list(pip_threshold = 0.5),
    debug = FALSE
  )
  config_b <- config_a
  config_b$parameters$pip_threshold <- 0.9

  hash_a <- cascade:::generate_config_hash(config_a)
  hash_b <- cascade:::generate_config_hash(config_b)
  expect_false(hash_a == hash_b)

  # Adding a threshold field that is NOT in the hash key list should NOT change hash
  config_c <- config_a
  config_c$thresholds <- list(lfsr_sig = 0.01)
  hash_c <- cascade:::generate_config_hash(config_c)
  expect_equal(hash_a, hash_c)
})

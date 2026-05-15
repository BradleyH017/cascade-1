# test-validation-genes.R
#
# Production gene validation: loads fixtures generated from raw production data,
# re-derives expected category via oracle, AND compares against C++ package output.
# Gated by CASCADE_TEST_LEVEL=slow.

skip_if(
  !identical(Sys.getenv("CASCADE_TEST_LEVEL", "fast"), "slow"),
  "Skipping slow validation test (set CASCADE_TEST_LEVEL=slow to run)"
)

fixture_file <- test_path("fixtures", "validation_genes.rds")
skip_if_not(file.exists(fixture_file), "Gene validation fixtures not generated")

fixtures <- readRDS(fixture_file)

test_that("all gene fixtures have required fields", {
  for (f in fixtures) {
    expect_true(!is.null(f$gene_id), info = "Missing gene_id")
    expect_true(!is.null(f$acat_qvalues), info = paste("Missing acat_qvalues for", f$gene_id))
    expect_true(!is.null(f$expected_category), info = paste("Missing expected_category for", f$gene_id))
    expect_true(f$expected_category %in% ORACLE_CATEGORIES,
      info = paste("Invalid category for", f$gene_id, ":", f$expected_category)
    )
  }
})

# Helper: build acat_matrix + LFSR lookup and call C++ categorize_features_from_acat
run_cpp_categorization <- function(fixture) {
  acat <- fixture$acat_qvalues
  wide <- data.table::dcast(acat, . ~ cell_type, value.var = "q_value")
  wide[, `.` := NULL]
  wide[, feature_id := fixture$gene_id]

  # Build LFSR lookup from fixture data (new long format!)
  lfsr_lookup <- NULL
  if (!is.null(fixture$lfsr_values) && nrow(fixture$lfsr_values) > 0) {
    lfsr_lookup <- data.table::data.table(
      feature_id = fixture$gene_id,
      cell_type = fixture$lfsr_values$cell_type,
      min_lfsr = fixture$lfsr_values$lfsr
    )
    # Drop NAs to match pipeline behavior
    lfsr_lookup <- lfsr_lookup[!is.na(min_lfsr)]
  }

  cpp_params <- cascade:::hierarchy_to_cpp_params(cascade::DEFAULT_CELL_HIERARCHY)
  result <- cascade:::categorize_features_from_acat(
    acat_matrix = wide,
    lfsr_lookup = lfsr_lookup,
    feature_type = "gene",
    l2_to_l1_mapping = cpp_params$l2_to_l1_mapping,
    lineage_groups = cpp_params$lineage_groups,
    subgroup_levels = cpp_params$subgroup_levels,
    bulk_cts = cpp_params$bulk_cts,
    other_cts = cpp_params$other_cts,
    sig_threshold = 0.05,
    lfsr_null_threshold = 0.5,
    specificity_categories = cpp_params$specificity_categories
  )
  result$primary_category[1]
}

for (i in seq_along(fixtures)) {
  fixture <- fixtures[[i]]

  test_that(paste0("gene ", fixture$gene_id, " → ", fixture$expected_category), {
    # Step 1: Re-derive from raw q-values via oracle
    sig_threshold <- 0.05
    sig_cts <- fixture$acat_qvalues[q_value < sig_threshold]$cell_type
    sig_l1 <- unique(oracle_map_to_l1(sig_cts))
    tested_l1 <- unique(oracle_map_to_l1(fixture$acat_qvalues$cell_type))

    lfsr_named <- NULL
    if (!is.null(fixture$lfsr_values) && nrow(fixture$lfsr_values) > 0) {
      lfsr_named <- setNames(fixture$lfsr_values$lfsr, fixture$lfsr_values$cell_type)
    }

    oracle_result <- oracle_categorize_specificity(
      sig_l1_cts = sig_l1,
      lfsr_values = lfsr_named,
      tested_l1_cts = tested_l1
    )

    # Step 2: Oracle must agree with stored expected
    expect_equal(oracle_result, fixture$expected_category,
      info = paste("Oracle re-derivation mismatch for", fixture$gene_id)
    )

    # Step 3: C++ package output (with LFSR) must agree with oracle (with LFSR)
    cpp_result <- run_cpp_categorization(fixture)
    expect_equal(cpp_result, oracle_result,
      info = paste(
        "C++ mismatch for", fixture$gene_id,
        "- C++:", cpp_result, "oracle:", oracle_result
      )
    )
  })
}

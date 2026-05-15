# test-integration-regression.R
#
# Integration and regression tests:
# - Output schema and cross-column invariant checks
# - Golden file regression
# - Output snapshot regression
#
# Gated by CASCADE_TEST_LEVEL: schema checks are medium, regression tests are slow.

library(data.table)

# ============================================================================
# SCHEMA AND INVARIANT CHECKS (medium tier)
# Applied to any available production/benchmark output
# ============================================================================

describe("Output schema and invariants", {
  skip_if(
    Sys.getenv("CASCADE_TEST_LEVEL", "fast") == "fast",
    "Skipping integration test (set CASCADE_TEST_LEVEL=medium or slow)"
  )

  # Set CASCADE_TEST_DATA_DIR (e.g. in ~/.Renviron) to a directory containing
  # variant_categorization.tsv.gz / gene_categorization.tsv.gz / peak_categorization.tsv.gz.
  test_data_dir <- Sys.getenv("CASCADE_TEST_DATA_DIR", "")
  skip_if(!nzchar(test_data_dir), "Set CASCADE_TEST_DATA_DIR to enable")
  skip_if_not(dir.exists(test_data_dir), "CASCADE_TEST_DATA_DIR not found")

  variant_file <- file.path(test_data_dir, "variant_categorization.tsv.gz")
  gene_file <- file.path(test_data_dir, "gene_categorization.tsv.gz")
  peak_file <- file.path(test_data_dir, "peak_categorization.tsv.gz")

  skip_if_not(file.exists(variant_file), "Variant output not found")

  v <- fread(variant_file)
  g <- fread(gene_file)
  p <- fread(peak_file)

  it("variant output has all expected columns", {
    expected_cols <- c(
      "variant_id", "qtl_mechanism_category", "cell_type_specificity",
      "gene_cell_type_specificity", "peak_cell_type_specificity",
      "qtl_pattern_number", "qtl_pattern", "best_cell_types",
      "significant_cts", "gene_affected_cell_types",
      "peak_affected_cell_types", "associated_genes",
      "associated_peaks", "cascade_peak_genes",
      "peak_overlap", "caqtl", "peak_gene_link", "eqtl",
      "link_mechanism"
    )
    for (col in expected_cols) {
      expect_true(col %in% names(v), info = paste("Missing column:", col))
    }
  })

  it("all mechanism categories are valid", {
    valid_mechs <- cascade::QTL_MECHANISMS
    invalid <- v[!qtl_mechanism_category %in% valid_mechs]
    expect_equal(nrow(invalid), 0,
      info = paste("Invalid mechanisms:", paste(unique(invalid$qtl_mechanism_category), collapse = ", "))
    )
  })

  it("all specificity categories are valid", {
    valid_specs <- cascade::DEFAULT_CELL_HIERARCHY$category_labels
    invalid <- v[!cell_type_specificity %in% valid_specs]
    expect_equal(nrow(invalid), 0,
      info = paste("Invalid specificities:", paste(unique(invalid$cell_type_specificity), collapse = ", "))
    )
  })

  it("link-dependent mechanism categories are present", {
    mechs <- unique(v$qtl_mechanism_category)
    # At least one cascade/link category should exist when both eQTL+caQTL data present
    link_mechs <- c("Local Cascade", "Positional Cascade", "Distal Cascade", "Only caQTL (With Link)")
    expect_true(any(link_mechs %in% mechs),
      info = paste("No link-dependent mechanisms found. Present:", paste(mechs, collapse = ", "))
    )
  })

  it("peak_gene_link is not 100% 'No link'", {
    pgl_vals <- unique(v$peak_gene_link)
    expect_true(length(pgl_vals) > 1 || pgl_vals[1] != "No link",
      info = "All variants have 'No link' — peak-gene link detection may be broken"
    )
  })

  it("eqtl field contains 'For linked gene' values", {
    expect_true("For linked gene" %in% v$eqtl,
      info = "No 'For linked gene' eQTL values — link detection may be broken"
    )
  })

  it("pattern numbers are all 1-22 (No molQTL filtered by PIP)", {
    expect_true(all(v$qtl_pattern_number >= 1 & v$qtl_pattern_number <= 22),
      info = paste("Patterns outside 1-22:", paste(unique(v$qtl_pattern_number[v$qtl_pattern_number > 22]), collapse = ", "))
    )
  })

  it("pattern → mechanism mapping is consistent", {
    for (i in seq_len(nrow(v))) {
      expected_mech <- oracle_pattern_to_mechanism(v$qtl_pattern_number[i])
      expect_equal(v$qtl_mechanism_category[i], expected_mech,
        info = paste(
          "Row", i, ": pattern", v$qtl_pattern_number[i],
          "maps to", expected_mech, "but got", v$qtl_mechanism_category[i]
        )
      )
    }
  })

  it("variant_id is never NA or empty", {
    expect_true(all(!is.na(v$variant_id) & v$variant_id != ""))
  })

  it("best_cell_types is never empty for QTL variants", {
    expect_true(all(!is.na(v$best_cell_types) & v$best_cell_types != ""))
  })

  it("peak_overlap is always TRUE or FALSE", {
    expect_true(all(v$peak_overlap %in% c(TRUE, FALSE)))
  })

  it("Local Cascade invariants hold", {
    fc <- v[qtl_mechanism_category == "Local Cascade"]
    if (nrow(fc) > 0) {
      expect_true(all(fc$peak_gene_link == "For caQTL peak"),
        info = "Local Cascade with wrong peak_gene_link"
      )
      expect_true(all(fc$eqtl == "For linked gene"),
        info = "Local Cascade with wrong eqtl"
      )
      expect_true(all(fc$associated_genes != "" & !is.na(fc$associated_genes)),
        info = "Local Cascade with empty associated_genes"
      )
      expect_true(all(fc$associated_peaks != "" & !is.na(fc$associated_peaks)),
        info = "Local Cascade with empty associated_peaks"
      )
    }
  })

  it("Positional Cascade invariants hold", {
    pc <- v[qtl_mechanism_category == "Positional Cascade"]
    if (nrow(pc) > 0) {
      # Positional Cascade: P2, P3 have overlap link; P4 has caQTL peak link
      expect_true(all(pc$peak_gene_link %in% c("For overlapping peak", "For caQTL peak")),
        info = "Positional Cascade with unexpected peak_gene_link"
      )
      expect_true(all(pc$eqtl == "For linked gene"),
        info = "Positional Cascade with wrong eqtl"
      )
      expect_true(all(pc$associated_genes != "" & !is.na(pc$associated_genes)),
        info = "Positional Cascade with empty associated_genes"
      )
    }
  })

  it("Distal Cascade invariants hold", {
    dc <- v[qtl_mechanism_category == "Distal Cascade"]
    if (nrow(dc) > 0) {
      expect_true(all(dc$peak_gene_link == "For caQTL peak"),
        info = "Distal Cascade with wrong peak_gene_link"
      )
      expect_true(all(dc$eqtl == "For linked gene"),
        info = "Distal Cascade with wrong eqtl"
      )
      expect_true(all(dc$associated_genes != "" & !is.na(dc$associated_genes)),
        info = "Distal Cascade with empty associated_genes"
      )
    }
  })

  it("Only caQTL (With Link) invariants hold", {
    cl <- v[qtl_mechanism_category == "Only caQTL (With Link)"]
    if (nrow(cl) > 0) {
      # P14 has 'For overlapping peak'; P13/15/16 have 'For caQTL peak'
      expect_true(all(cl$peak_gene_link %in% c("For caQTL peak", "For overlapping peak")),
        info = "caQTL With Link with unexpected peak_gene_link"
      )
      expect_true(all(cl$eqtl == "No detected eQTL"),
        info = "caQTL With Link with detected eQTL"
      )
    }
  })

  it("gene output has expected columns", {
    expect_true("gene_id" %in% names(g), info = "gene output missing gene_id")
    expect_true("cell_type_specificity" %in% names(g), info = "gene output missing cell_type_specificity")
    expect_true("n_significant_cts" %in% names(g), info = "gene output missing n_significant_cts")
    expect_true("top_variant" %in% names(g), info = "gene output missing top_variant")
  })

  it("peak output has expected columns", {
    expect_true("peak_id" %in% names(p), info = "peak output missing peak_id")
    expect_true("cell_type_specificity" %in% names(p), info = "peak output missing cell_type_specificity")
    expect_true("n_significant_cts" %in% names(p), info = "peak output missing n_significant_cts")
  })
})

# ============================================================================
# GOLDEN FILE REGRESSION (slow tier)
# ============================================================================

describe("Golden file regression", {
  skip_if(
    !identical(Sys.getenv("CASCADE_TEST_LEVEL", "fast"), "slow"),
    "Skipping slow regression test"
  )

  golden_file <- test_path("fixtures", "golden_variants.csv")
  skip_if_not(file.exists(golden_file), "Golden file not generated")

  test_data_dir <- Sys.getenv("CASCADE_TEST_DATA_DIR", "")
  skip_if(!nzchar(test_data_dir), "Set CASCADE_TEST_DATA_DIR to enable")
  variant_file <- file.path(test_data_dir, "variant_categorization.tsv.gz")
  skip_if_not(file.exists(variant_file), "variant_categorization.tsv.gz not found in CASCADE_TEST_DATA_DIR")

  golden <- fread(golden_file)
  current <- fread(variant_file)

  it("golden file variants exist in current output", {
    missing <- setdiff(golden$variant_id, current$variant_id)
    expect_equal(length(missing), 0,
      info = paste("Missing variants:", paste(head(missing, 5), collapse = ", "))
    )
  })

  it("golden file rows match current output exactly", {
    stable_cols <- names(golden)
    current_subset <- current[variant_id %in% golden$variant_id][order(variant_id)]
    current_subset <- current_subset[, ..stable_cols]

    for (col in stable_cols) {
      mismatches <- golden[[col]] != current_subset[[col]]
      mismatches[is.na(golden[[col]]) & is.na(current_subset[[col]])] <- FALSE
      n_mismatch <- sum(mismatches, na.rm = TRUE)
      expect_equal(n_mismatch, 0,
        info = paste("Column", col, "has", n_mismatch, "mismatches")
      )
    }
  })
})

# ============================================================================
# OUTPUT SNAPSHOT REGRESSION (slow tier)
# ============================================================================

describe("Output snapshot regression", {
  skip_if(
    !identical(Sys.getenv("CASCADE_TEST_LEVEL", "fast"), "slow"),
    "Skipping slow regression test"
  )

  test_data_dir <- Sys.getenv("CASCADE_TEST_DATA_DIR", "")
  skip_if(!nzchar(test_data_dir), "Set CASCADE_TEST_DATA_DIR to enable")
  skip_if_not(dir.exists(test_data_dir), "CASCADE_TEST_DATA_DIR not found")

  v <- fread(file.path(test_data_dir, "variant_categorization.tsv.gz"))

  it("variant count matches snapshot", {
    expect_equal(nrow(v), 3219)
  })

  it("mechanism distribution matches snapshot", {
    mech_dist <- v[, .N, by = qtl_mechanism_category]
    setkey(mech_dist, qtl_mechanism_category)
    expect_equal(mech_dist["Only caQTL (With Link)", N], 1684)
    expect_equal(mech_dist["Only caQTL (No Link)", N], 1145)
    expect_equal(mech_dist["Only eQTL", N], 217)
    expect_equal(mech_dist["Local Cascade", N], 37)
    expect_equal(mech_dist["Positional Cascade", N], 38)
    expect_equal(mech_dist["Distal Cascade", N], 28)
    expect_equal(mech_dist["caQTL + eQTL (No Link)", N], 70)
  })

  it("specificity distribution matches snapshot", {
    spec_dist <- v[, .N, by = cell_type_specificity]
    setkey(spec_dist, cell_type_specificity)
    expect_equal(spec_dist["Single cell-type", N], 1751)
    expect_equal(spec_dist["Likely shared but underpowered", N], 535)
    expect_equal(spec_dist["Cross-lineage shared", N], 409)
    expect_equal(spec_dist["Lineage-specific", N], 305)
    expect_equal(spec_dist["T-cell-specific", N], 219)
  })
})

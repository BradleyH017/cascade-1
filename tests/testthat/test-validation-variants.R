# test-validation-variants.R
#
# Production variant validation: loads fixtures with per-CT SuSiE data,
# verifies production output matches expected mechanism, specificity,
# and cross-column invariants.
# Gated by CASCADE_TEST_LEVEL=slow.

skip_if(
  !identical(Sys.getenv("CASCADE_TEST_LEVEL", "fast"), "slow"),
  "Skipping slow validation test (set CASCADE_TEST_LEVEL=slow to run)"
)

fixture_file <- test_path("fixtures", "validation_variants.rds")
skip_if_not(file.exists(fixture_file), "Variant validation fixtures not generated")

fixtures <- readRDS(fixture_file)

test_that("all variant fixtures have required fields", {
  for (f in fixtures) {
    expect_true(!is.null(f$variant_id), info = "Missing variant_id")
    expect_true(!is.null(f$per_celltype), info = paste("Missing per_celltype for", f$variant_id))
    expect_true(!is.null(f$production_output), info = paste("Missing production_output for", f$variant_id))
    expect_true(f$production_output$qtl_mechanism_category %in% cascade::QTL_MECHANISMS,
      info = paste("Invalid mechanism for", f$variant_id)
    )
  }
})

for (i in seq_along(fixtures)) {
  fixture <- fixtures[[i]]

  test_that(paste0("variant ", fixture$variant_id, " → ", fixture$production_output$qtl_mechanism_category), {
    prod <- fixture$production_output
    pip_threshold <- 0.01

    # 1. Pattern → mechanism consistency (via oracle)
    expected_mechanism <- oracle_pattern_to_mechanism(prod$qtl_pattern_number)
    expect_equal(prod$qtl_mechanism_category, expected_mechanism,
      info = paste(
        "Pattern", prod$qtl_pattern_number,
        "should map to", expected_mechanism
      )
    )

    # 2. Specificity is valid
    expect_true(prod$cell_type_specificity %in% cascade::DEFAULT_CELL_HIERARCHY$category_labels,
      info = paste("Invalid specificity:", prod$cell_type_specificity)
    )

    # 3. Per-CT data consistency
    cts_with_eqtl <- names(Filter(function(ct) ct$eqtl_pip > pip_threshold, fixture$per_celltype))
    cts_with_caqtl <- names(Filter(function(ct) ct$caqtl_pip > pip_threshold, fixture$per_celltype))

    # 4. Mechanism ↔ field cross-column invariants
    mech <- prod$qtl_mechanism_category
    patt <- prod$qtl_pattern_number

    # Local Cascade (P1): caQTL for overlapping peak + eQTL for linked gene + overlap link
    if (mech == "Local Cascade") {
      expect_true(length(cts_with_eqtl) > 0,
        info = paste(fixture$variant_id, "Local Cascade needs eQTL CTs")
      )
      expect_true(length(cts_with_caqtl) > 0,
        info = paste(fixture$variant_id, "Local Cascade needs caQTL CTs")
      )
      expect_equal(prod$peak_gene_link, "For caQTL peak",
        info = paste(fixture$variant_id, "Local Cascade link must be 'For caQTL peak'")
      )
      expect_equal(prod$eqtl, "For linked gene",
        info = paste(fixture$variant_id, "Local Cascade eqtl must be 'For linked gene'")
      )
    }

    # Positional Cascade (P2-4): link present + eQTL for linked gene
    # link can be 'For overlapping peak' (P2, P3) or 'For caQTL peak' (P4)
    if (mech == "Positional Cascade") {
      expect_true(prod$peak_gene_link %in% c("For overlapping peak", "For caQTL peak"),
        info = paste(fixture$variant_id, "Positional Cascade link =", prod$peak_gene_link)
      )
      expect_equal(prod$eqtl, "For linked gene",
        info = paste(fixture$variant_id, "Positional Cascade eqtl must be 'For linked gene'")
      )
    }

    # Distal Cascade (P5): no overlap + non-overlap caQTL + non-overlap link + eQTL linked
    if (mech == "Distal Cascade") {
      expect_true(length(cts_with_eqtl) > 0,
        info = paste(fixture$variant_id, "Distal Cascade needs eQTL CTs")
      )
      expect_true(length(cts_with_caqtl) > 0,
        info = paste(fixture$variant_id, "Distal Cascade needs caQTL CTs")
      )
      expect_equal(prod$peak_gene_link, "For caQTL peak",
        info = paste(fixture$variant_id, "Distal Cascade link must be 'For caQTL peak'")
      )
      expect_equal(prod$eqtl, "For linked gene",
        info = paste(fixture$variant_id, "Distal Cascade eqtl must be 'For linked gene'")
      )
    }

    # Only caQTL (With Link) (P13-16): link present, no eQTL, must have caQTL
    # P13/15/16 have 'For caQTL peak'; P14 has 'For overlapping peak'
    if (mech == "Only caQTL (With Link)") {
      expect_true(length(cts_with_caqtl) > 0,
        info = paste(fixture$variant_id, "caQTL With Link must have caQTL CTs")
      )
      expect_true(prod$peak_gene_link %in% c("For caQTL peak", "For overlapping peak"),
        info = paste(fixture$variant_id, "caQTL With Link link =", prod$peak_gene_link)
      )
      expect_equal(prod$eqtl, "No detected eQTL",
        info = paste(fixture$variant_id, "caQTL With Link must have no eQTL")
      )
    }

    # Only caQTL (No Link) (P17-19): no link, no eQTL, must have caQTL
    if (mech == "Only caQTL (No Link)") {
      expect_true(length(cts_with_caqtl) > 0,
        info = paste(fixture$variant_id, "caQTL No Link must have caQTL CTs")
      )
      expect_equal(prod$peak_gene_link, "No link",
        info = paste(fixture$variant_id, "caQTL No Link: link =", prod$peak_gene_link)
      )
      expect_equal(prod$eqtl, "No detected eQTL",
        info = paste(fixture$variant_id, "caQTL No Link must have no eQTL")
      )
    }

    # Only eQTL (P20-22): must have eQTL evidence
    if (mech == "Only eQTL") {
      expect_true(length(cts_with_eqtl) > 0,
        info = paste(fixture$variant_id, "Only eQTL needs eQTL CTs")
      )
      expect_true(prod$eqtl %in% c("eQTL", "For non-linked gene"),
        info = paste(fixture$variant_id, "Only eQTL: eqtl =", prod$eqtl)
      )
    }

    # caQTL + eQTL (No Link) (P6-12): BOTH caQTL and eQTL present but not linked
    if (mech == "caQTL + eQTL (No Link)") {
      expect_true(length(cts_with_eqtl) > 0,
        info = paste(fixture$variant_id, "caQTL+eQTL must have eQTL CTs")
      )
      expect_true(length(cts_with_caqtl) > 0,
        info = paste(fixture$variant_id, "caQTL+eQTL must have caQTL CTs")
      )
      expect_true(prod$eqtl %in% c("eQTL", "For non-linked gene"),
        info = paste(fixture$variant_id, "caQTL+eQTL: eqtl =", prod$eqtl)
      )
    }
  })
}

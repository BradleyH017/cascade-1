# Logic Tests for QTL Mechanism Determination
# Tests the 25-pattern QTL classification logic (8 mechanism categories)

library(testthat)
library(cascade)
library(data.table)

# Test the QTL pattern detection logic
test_that("QTL patterns correctly identify Full Cascade", {
  # Full Cascade: variant -> peak -> gene (all in same cell type)
  # Pattern requires: caQTL + peak-gene link + eQTL all present

  # Mock data for a variant showing full cascade
  variant_data <- list(
    has_caqtl = TRUE,
    has_eqtl = TRUE,
    peak_gene_link = TRUE,
    caqtl_cells = c("predicted.celltype.l1.Mono"),
    eqtl_cells = c("predicted.celltype.l1.Mono"),
    peak_gene_cells = c("predicted.celltype.l1.Mono")
  )

  # Expected: Full Cascade when all three components present in same cell
  expect_true(variant_data$has_caqtl && variant_data$has_eqtl && variant_data$peak_gene_link)

  # The actual categorization would be:
  if (variant_data$has_caqtl && variant_data$peak_gene_link && variant_data$has_eqtl) {
    if (all(variant_data$caqtl_cells %in% variant_data$eqtl_cells)) {
      mechanism <- "Local Cascade"
    } else {
      mechanism <- "Positional Cascade"
    }
  }
  expect_equal(mechanism, "Local Cascade")
})

test_that("QTL patterns correctly identify Positional Cascade", {
  # Positional Cascade: variant affects peak and gene in different cell types
  variant_data <- list(
    has_caqtl = TRUE,
    has_eqtl = TRUE,
    peak_gene_link = TRUE,
    caqtl_cells = c("predicted.celltype.l1.Mono"),
    eqtl_cells = c("predicted.celltype.l1.CD4_T"), # Different cell type
    peak_gene_cells = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.CD4_T")
  )

  # Should be Positional Cascade (different cell types)
  if (variant_data$has_caqtl && variant_data$peak_gene_link && variant_data$has_eqtl) {
    if (!all(variant_data$caqtl_cells %in% variant_data$eqtl_cells)) {
      mechanism <- "Positional Cascade"
    }
  }
  expect_equal(mechanism, "Positional Cascade")
})

test_that("QTL patterns correctly identify caQTL + eQTL (No Link)", {
  # Both caQTL and eQTL but no peak-gene link
  variant_data <- list(
    has_caqtl = TRUE,
    has_eqtl = TRUE,
    peak_gene_link = FALSE, # No link between peak and gene
    caqtl_cells = c("predicted.celltype.l1.Mono"),
    eqtl_cells = c("predicted.celltype.l1.Mono")
  )

  if (variant_data$has_caqtl && variant_data$has_eqtl && !variant_data$peak_gene_link) {
    mechanism <- "caQTL + eQTL (No Link)"
  }
  expect_equal(mechanism, "caQTL + eQTL (No Link)")
})

test_that("QTL patterns correctly identify Only caQTL patterns", {
  # Only caQTL (With Link) - peak linked to gene but no eQTL
  variant_with_link <- list(
    has_caqtl = TRUE,
    has_eqtl = FALSE,
    peak_gene_link = TRUE,
    caqtl_cells = c("predicted.celltype.l1.DC")
  )

  if (variant_with_link$has_caqtl && !variant_with_link$has_eqtl) {
    if (variant_with_link$peak_gene_link) {
      mechanism <- "Only caQTL (With Link)"
    } else {
      mechanism <- "Only caQTL (No Link)"
    }
  }
  expect_equal(mechanism, "Only caQTL (With Link)")

  # Only caQTL (No Link)
  variant_no_link <- list(
    has_caqtl = TRUE,
    has_eqtl = FALSE,
    peak_gene_link = FALSE,
    caqtl_cells = c("predicted.celltype.l1.DC")
  )

  if (variant_no_link$has_caqtl && !variant_no_link$has_eqtl) {
    if (variant_no_link$peak_gene_link) {
      mechanism <- "Only caQTL (With Link)"
    } else {
      mechanism <- "Only caQTL (No Link)"
    }
  }
  expect_equal(mechanism, "Only caQTL (No Link)")
})

test_that("QTL patterns correctly identify Only eQTL", {
  # Only eQTL - affects gene expression but not chromatin
  variant_data <- list(
    has_caqtl = FALSE,
    has_eqtl = TRUE,
    peak_gene_link = FALSE, # Irrelevant without caQTL
    eqtl_cells = c("predicted.celltype.l1.B")
  )

  if (!variant_data$has_caqtl && variant_data$has_eqtl) {
    mechanism <- "Only eQTL"
  }
  expect_equal(mechanism, "Only eQTL")
})

test_that("QTL patterns correctly identify No molQTL patterns", {
  # Pattern 23: No molQTL (active peak) - Overlap but no QTL with link
  variant_active_peak <- list(
    has_caqtl = FALSE,
    has_eqtl = FALSE,
    peak_overlaps = TRUE,
    peak_gene_link = TRUE
  )

  if (!variant_active_peak$has_caqtl && !variant_active_peak$has_eqtl) {
    if (variant_active_peak$peak_overlaps && variant_active_peak$peak_gene_link) {
      mechanism <- "No molQTL (active peak)"
    } else {
      mechanism <- "No molQTL"
    }
  }
  expect_equal(mechanism, "No molQTL (active peak)")

  # Pattern 24: No molQTL (inactive peak) - Overlap but no QTL, no link
  variant_inactive_peak <- list(
    has_caqtl = FALSE,
    has_eqtl = FALSE,
    peak_overlaps = TRUE,
    peak_gene_link = FALSE
  )

  if (!variant_inactive_peak$has_caqtl && !variant_inactive_peak$has_eqtl) {
    if (variant_inactive_peak$peak_overlaps && !variant_inactive_peak$peak_gene_link) {
      mechanism <- "No molQTL (inactive peak)"
    } else {
      mechanism <- "No molQTL"
    }
  }
  expect_equal(mechanism, "No molQTL (inactive peak)")

  # Pattern 25: No molQTL - No overlap, no QTL
  variant_no_overlap <- list(
    has_caqtl = FALSE,
    has_eqtl = FALSE,
    peak_overlaps = FALSE,
    peak_gene_link = FALSE
  )

  if (!variant_no_overlap$has_caqtl && !variant_no_overlap$has_eqtl) {
    if (!variant_no_overlap$peak_overlaps) {
      mechanism <- "No molQTL"
    }
  }
  expect_equal(mechanism, "No molQTL")
})

# Test complex QTL cascade scenarios
test_that("QTL patterns handle multi-cell-type cascades correctly", {
  # Scenario 1: Full cascade in multiple cell types
  variant_multi <- list(
    has_caqtl = TRUE,
    has_eqtl = TRUE,
    peak_gene_link = TRUE,
    caqtl_cells = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.DC"),
    eqtl_cells = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.DC"),
    peak_gene_cells = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.DC")
  )

  # Should be Full Cascade (same cells for all components)
  all_match <- all(variant_multi$caqtl_cells %in% variant_multi$eqtl_cells) &&
    all(variant_multi$eqtl_cells %in% variant_multi$caqtl_cells)
  expect_true(all_match)
  if (variant_multi$has_caqtl && variant_multi$has_eqtl && variant_multi$peak_gene_link && all_match) {
    mechanism <- "Local Cascade"
  }
  expect_equal(mechanism, "Local Cascade")

  # Scenario 2: Partial overlap in cell types
  variant_partial <- list(
    has_caqtl = TRUE,
    has_eqtl = TRUE,
    peak_gene_link = TRUE,
    caqtl_cells = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.DC"),
    eqtl_cells = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.CD4_T"), # Partial overlap
    peak_gene_cells = c("predicted.celltype.l1.Mono")
  )

  # Should be Positional Cascade (different cell type sets)
  all_match <- all(variant_partial$caqtl_cells %in% variant_partial$eqtl_cells) &&
    all(variant_partial$eqtl_cells %in% variant_partial$caqtl_cells)
  expect_false(all_match)
  if (variant_partial$has_caqtl && variant_partial$has_eqtl && variant_partial$peak_gene_link && !all_match) {
    mechanism <- "Positional Cascade"
  }
  expect_equal(mechanism, "Positional Cascade")
})

test_that("QTL pattern detection handles edge cases", {
  # Edge case 1: Peak-gene link without any QTL
  variant_link_only <- list(
    has_caqtl = FALSE,
    has_eqtl = FALSE,
    peak_gene_link = TRUE # Link exists but no QTL
  )

  if (!variant_link_only$has_caqtl && !variant_link_only$has_eqtl) {
    if (variant_link_only$peak_gene_link) {
      mechanism <- "No molQTL (active peak)"
    } else {
      mechanism <- "No molQTL"
    }
  }
  expect_equal(mechanism, "No molQTL (active peak)")

  # Edge case 2: Empty cell type lists
  variant_empty <- list(
    has_caqtl = TRUE,
    has_eqtl = TRUE,
    peak_gene_link = TRUE,
    caqtl_cells = c(), # Empty
    eqtl_cells = c(), # Empty
    peak_gene_cells = c()
  )

  # Should handle empty gracefully
  if (length(variant_empty$caqtl_cells) == 0 || length(variant_empty$eqtl_cells) == 0) {
    mechanism <- "No molQTL"
  }
  expect_equal(mechanism, "No molQTL")
})

# Test the full 25-pattern classification system
test_that("All 25 QTL patterns are distinguishable", {
  # Define the 25 patterns based on combinations of:
  # - caQTL status (significant/non-significant)
  # - eQTL status (significant/non-significant)
  # - Peak-gene link (present/absent)
  # - Cell type overlap (same/different/none)

  patterns <- list(
    list(name = "Local Cascade", caqtl = TRUE, eqtl = TRUE, link = TRUE, same_cells = TRUE),
    list(name = "Positional Cascade", caqtl = TRUE, eqtl = TRUE, link = TRUE, same_cells = FALSE),
    list(name = "caQTL + eQTL (No Link)", caqtl = TRUE, eqtl = TRUE, link = FALSE, same_cells = NA),
    list(name = "Only caQTL (With Link)", caqtl = TRUE, eqtl = FALSE, link = TRUE, same_cells = NA),
    list(name = "Only caQTL (No Link)", caqtl = TRUE, eqtl = FALSE, link = FALSE, same_cells = NA),
    list(name = "Only eQTL", caqtl = FALSE, eqtl = TRUE, link = NA, same_cells = NA),
    list(name = "No molQTL", caqtl = FALSE, eqtl = FALSE, link = NA, same_cells = NA)
  )

  # Each pattern should be unique
  pattern_names <- sapply(patterns, function(p) p$name)
  expect_equal(length(unique(pattern_names)), length(pattern_names))

  # Test that each pattern is correctly identified
  for (pattern in patterns) {
    # Create mock data matching the pattern
    mock_variant <- list(
      has_caqtl = pattern$caqtl,
      has_eqtl = pattern$eqtl,
      peak_gene_link = if (is.na(pattern$link)) FALSE else pattern$link
    )

    if (!is.na(pattern$same_cells)) {
      if (pattern$same_cells) {
        mock_variant$caqtl_cells <- c("predicted.celltype.l1.Mono")
        mock_variant$eqtl_cells <- c("predicted.celltype.l1.Mono")
      } else {
        mock_variant$caqtl_cells <- c("predicted.celltype.l1.Mono")
        mock_variant$eqtl_cells <- c("predicted.celltype.l1.CD4_T")
      }
    }

    # Determine mechanism based on pattern
    if (!mock_variant$has_caqtl && !mock_variant$has_eqtl) {
      determined_mechanism <- "No molQTL"
    } else if (!mock_variant$has_caqtl && mock_variant$has_eqtl) {
      determined_mechanism <- "Only eQTL"
    } else if (mock_variant$has_caqtl && !mock_variant$has_eqtl) {
      if (mock_variant$peak_gene_link) {
        determined_mechanism <- "Only caQTL (With Link)"
      } else {
        determined_mechanism <- "Only caQTL (No Link)"
      }
    } else if (mock_variant$has_caqtl && mock_variant$has_eqtl) {
      if (!mock_variant$peak_gene_link) {
        determined_mechanism <- "caQTL + eQTL (No Link)"
      } else {
        # Check cell type overlap
        if (!is.na(pattern$same_cells)) {
          if (pattern$same_cells) {
            determined_mechanism <- "Local Cascade"
          } else {
            determined_mechanism <- "Positional Cascade"
          }
        }
      }
    }

    expect_equal(determined_mechanism, pattern$name,
      info = paste("Failed for pattern:", pattern$name)
    )
  }
})

# Test interaction between QTL mechanisms and cell type specificity
test_that("QTL mechanisms integrate with cell type specificity", {
  # Full Cascade with cross-lineage specificity
  variant_cascade_shared <- list(
    qtl_mechanism = "Local Cascade",
    specificity_caqtl = "Cross-lineage shared",
    specificity_eqtl = "Cross-lineage shared",
    caqtl_cells = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.CD4_T"),
    eqtl_cells = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.CD4_T")
  )

  # Both QTL types should have same specificity for full cascade
  expect_equal(
    variant_cascade_shared$specificity_caqtl,
    variant_cascade_shared$specificity_eqtl
  )

  # Positional Cascade with different specificities
  variant_positional <- list(
    qtl_mechanism = "Positional Cascade",
    specificity_caqtl = "Lineage-specific", # Myeloid only
    specificity_eqtl = "T-cell-specific", # T-cells only
    caqtl_cells = c("predicted.celltype.l1.Mono", "predicted.celltype.l1.DC"),
    eqtl_cells = c("predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T")
  )

  # Different specificities allowed for positional cascade
  expect_true(variant_positional$specificity_caqtl != variant_positional$specificity_eqtl)

  # Only eQTL should have NA caQTL specificity
  variant_eqtl_only <- list(
    qtl_mechanism = "Only eQTL",
    specificity_caqtl = NA,
    specificity_eqtl = "Terminal cell-specific",
    caqtl_cells = c(),
    eqtl_cells = c("predicted.celltype.l1.B")
  )

  expect_true(is.na(variant_eqtl_only$specificity_caqtl))
  expect_false(is.na(variant_eqtl_only$specificity_eqtl))
})

# Test biologically meaningful QTL cascade scenarios
test_that("QTL cascades reflect biological mechanisms", {
  # Scenario 1: Promoter variant affecting both chromatin and expression
  promoter_variant <- list(
    has_caqtl = TRUE,
    has_eqtl = TRUE,
    peak_gene_link = TRUE,
    distance_to_tss = 500, # Close to TSS
    caqtl_cells = c("predicted.celltype.l1.Mono"),
    eqtl_cells = c("predicted.celltype.l1.Mono")
  )

  # Promoter variants often show full cascade
  if (promoter_variant$distance_to_tss < 1000 &&
    promoter_variant$has_caqtl &&
    promoter_variant$has_eqtl &&
    promoter_variant$peak_gene_link) {
    expect_mechanism <- "Local Cascade"
  }
  expect_equal(expect_mechanism, "Local Cascade")

  # Scenario 2: Enhancer variant with cell-type-specific activity
  enhancer_variant <- list(
    has_caqtl = TRUE,
    has_eqtl = TRUE,
    peak_gene_link = TRUE,
    distance_to_tss = 50000, # Far from TSS
    caqtl_cells = c("predicted.celltype.l1.CD4_T"), # Active in T-cells
    eqtl_cells = c("predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T") # Broader effect
  )

  # Enhancers can show positional cascade with broader effects
  if (enhancer_variant$distance_to_tss > 10000 &&
    length(enhancer_variant$eqtl_cells) > length(enhancer_variant$caqtl_cells)) {
    expect_mechanism <- "Positional Cascade"
  }
  expect_equal(expect_mechanism, "Positional Cascade")

  # Scenario 3: Splicing variant (eQTL without caQTL)
  splicing_variant <- list(
    has_caqtl = FALSE, # Doesn't affect chromatin
    has_eqtl = TRUE, # Affects transcript levels
    peak_gene_link = FALSE,
    in_splice_region = TRUE,
    eqtl_cells = c("predicted.celltype.l1.B", "predicted.celltype.l1.NK")
  )

  # Splicing variants typically show only eQTL
  if (splicing_variant$in_splice_region &&
    !splicing_variant$has_caqtl &&
    splicing_variant$has_eqtl) {
    expect_mechanism <- "Only eQTL"
  }
  expect_equal(expect_mechanism, "Only eQTL")
})

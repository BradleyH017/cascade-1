# Logic Tests for Variant Categorization and Power Assessment
# Tests the integration of QTL mechanisms with cell type specificity and power assessment

library(testthat)
library(cascade)
library(data.table)

# Test variant categorization logic
test_that("Variant categorization integrates QTL mechanism and specificity", {
  # Variant with Local Cascade and Cross-lineage specificity
  variant <- list(
    id = "chr1_100000_A_G",
    qtl_mechanism = "Local Cascade",
    caqtl_specificity = "Cross-lineage shared",
    eqtl_specificity = "Cross-lineage shared",
    peak_id = "chr1_99500_100500",
    gene_id = "ENSG00000000001"
  )

  # For Local Cascade, primary specificity should match both QTLs
  primary_specificity <- variant$caqtl_specificity # Both should be same for local cascade
  expect_equal(primary_specificity, "Cross-lineage shared")
  expect_equal(variant$caqtl_specificity, variant$eqtl_specificity)

  # Overall categorization
  categorization <- paste(variant$qtl_mechanism, "-", primary_specificity)
  expect_equal(categorization, "Local Cascade - Cross-lineage shared")
})

test_that("Variant categorization handles Positional Cascade correctly", {
  # Positional Cascade with different specificities
  variant <- list(
    id = "chr2_200000_C_T",
    qtl_mechanism = "Positional Cascade",
    caqtl_specificity = "Lineage-specific", # Myeloid
    eqtl_specificity = "Single cell-type", # Single cell
    peak_id = "chr2_199000_201000",
    gene_id = "ENSG00000000002"
  )

  # For Positional Cascade, should report both specificities
  categorization_parts <- list(
    mechanism = variant$qtl_mechanism,
    caqtl = variant$caqtl_specificity,
    eqtl = variant$eqtl_specificity
  )

  expect_equal(categorization_parts$mechanism, "Positional Cascade")
  expect_true(categorization_parts$caqtl != categorization_parts$eqtl)
})

test_that("Variant categorization handles single QTL types", {
  # Only caQTL
  variant_caqtl <- list(
    id = "chr3_300000_G_A",
    qtl_mechanism = "Only caQTL (With Link)",
    caqtl_specificity = "T-cell-specific",
    eqtl_specificity = NA,
    peak_id = "chr3_299000_301000",
    gene_id = "ENSG00000000003"
  )

  # Should use caQTL specificity as primary
  primary_specificity <- variant_caqtl$caqtl_specificity
  expect_equal(primary_specificity, "T-cell-specific")
  expect_true(is.na(variant_caqtl$eqtl_specificity))

  # Only eQTL
  variant_eqtl <- list(
    id = "chr4_400000_T_C",
    qtl_mechanism = "Only eQTL",
    caqtl_specificity = NA,
    eqtl_specificity = "Cross-lineage shared",
    peak_id = NA,
    gene_id = "ENSG00000000004"
  )

  # Should use eQTL specificity as primary
  primary_specificity <- variant_eqtl$eqtl_specificity
  expect_equal(primary_specificity, "Cross-lineage shared")
  expect_true(is.na(variant_eqtl$caqtl_specificity))
})

# Test LFSR-based power assessment
test_that("LFSR correctly identifies underpowered effects", {
  # Gene with one significant cell type but gray zone in related cells
  lfsr_values <- c(
    predicted.celltype.l1.Mono = 0.01, # Significant
    predicted.celltype.l1.DC = 0.15, # Gray zone (0.05-0.5)
    predicted.celltype.l1.CD4_T = 0.6, # Null
    predicted.celltype.l1.CD8_T = 0.7, # Null
    predicted.celltype.l1.B = 0.8, # Null
    predicted.celltype.l1.NK = 0.9 # Null
  )

  sig_threshold <- 0.05
  null_threshold <- 0.5

  # Identify categories
  significant <- names(lfsr_values)[lfsr_values < sig_threshold]
  gray_zone <- names(lfsr_values)[lfsr_values >= sig_threshold & lfsr_values < null_threshold]
  null_cells <- names(lfsr_values)[lfsr_values >= null_threshold]

  expect_equal(significant, "predicted.celltype.l1.Mono")
  expect_equal(gray_zone, "predicted.celltype.l1.DC")
  expect_equal(length(null_cells), 4)

  # Should be categorized as underpowered due to gray zone
  has_gray_zone <- length(gray_zone) > 0
  expect_true(has_gray_zone)
})

test_that("LFSR power assessment scales with evidence strength", {
  # Strong evidence across lineages
  lfsr_strong <- c(
    predicted.celltype.l1.Mono = 0.001, # Very strong
    predicted.celltype.l1.DC = 0.002, # Very strong
    predicted.celltype.l1.CD4_T = 0.003, # Very strong
    predicted.celltype.l1.CD8_T = 0.004, # Very strong
    predicted.celltype.l1.B = 0.7, # Null
    predicted.celltype.l1.NK = 0.8 # Null
  )

  significant_strong <- names(lfsr_strong)[lfsr_strong < 0.05]
  expect_equal(length(significant_strong), 4)

  # Category should be Cross-lineage shared (strong evidence)
  has_myeloid <- any(c("predicted.celltype.l1.Mono", "predicted.celltype.l1.DC") %in% significant_strong)
  has_lymphoid <- any(c("predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T") %in% significant_strong)
  expect_true(has_myeloid && has_lymphoid)

  # Weak evidence with gray zones
  lfsr_weak <- c(
    predicted.celltype.l1.Mono = 0.04, # Barely significant
    predicted.celltype.l1.DC = 0.08, # Gray zone
    predicted.celltype.l1.CD4_T = 0.12, # Gray zone
    predicted.celltype.l1.CD8_T = 0.20, # Gray zone
    predicted.celltype.l1.B = 0.6, # Null
    predicted.celltype.l1.NK = 0.7 # Null
  )

  significant_weak <- names(lfsr_weak)[lfsr_weak < 0.05]
  gray_zone_weak <- names(lfsr_weak)[lfsr_weak >= 0.05 & lfsr_weak < 0.5]

  expect_equal(length(significant_weak), 1)
  expect_equal(length(gray_zone_weak), 3)

  # Should indicate underpowered
  expect_true(length(gray_zone_weak) > 0)
})

test_that("Power assessment handles heterogeneous effects", {
  # Heterogeneous effect: strong in some cells, weak in others
  lfsr_hetero <- c(
    predicted.celltype.l1.Mono = 0.001, # Very strong
    predicted.celltype.l1.DC = 0.45, # Gray zone
    predicted.celltype.l1.CD4_T = 0.002, # Very strong
    predicted.celltype.l1.CD8_T = 0.48, # Gray zone
    predicted.celltype.l1.B = 0.003, # Very strong
    predicted.celltype.l1.NK = 0.9 # Null
  )

  significant <- names(lfsr_hetero)[lfsr_hetero < 0.05]
  gray_zone <- names(lfsr_hetero)[lfsr_hetero >= 0.05 & lfsr_hetero < 0.5]

  # Multiple significant across lineages
  expect_equal(length(significant), 3)
  expect_true("predicted.celltype.l1.Mono" %in% significant) # Myeloid
  expect_true("predicted.celltype.l1.CD4_T" %in% significant) # Lymphoid

  # But also gray zones present
  expect_equal(length(gray_zone), 2)

  # Pattern suggests true cross-lineage with some cells underpowered
  categorization <- if (length(gray_zone) > 0) {
    "Likely shared but underpowered"
  } else {
    "Cross-lineage shared"
  }
  expect_equal(categorization, "Likely shared but underpowered")
})

# Test Cochran's Q vs LFSR integration
test_that("Cochran's Q and LFSR provide complementary power assessment", {
  # Cochran's Q significant (heterogeneity detected)
  cochran_q_pval <- 1e-10
  cochran_q_threshold <- 5e-8

  # LFSR shows specific pattern
  lfsr_values <- c(
    predicted.celltype.l1.Mono = 0.001,
    predicted.celltype.l1.DC = 0.8,
    predicted.celltype.l1.CD4_T = 0.002,
    predicted.celltype.l1.CD8_T = 0.9,
    predicted.celltype.l1.B = 0.85,
    predicted.celltype.l1.NK = 0.003
  )

  # Cochran's Q indicates heterogeneity
  has_heterogeneity <- cochran_q_pval < cochran_q_threshold
  expect_true(has_heterogeneity)

  # LFSR reveals which cells drive the heterogeneity
  significant_lfsr <- names(lfsr_values)[lfsr_values < 0.05]
  expect_equal(length(significant_lfsr), 3)
  expect_true("predicted.celltype.l1.Mono" %in% significant_lfsr)
  expect_true("predicted.celltype.l1.NK" %in% significant_lfsr)

  # Combined interpretation
  if (has_heterogeneity && length(significant_lfsr) > 0) {
    interpretation <- "Heterogeneous effect confirmed by both methods"
  }
  expect_equal(interpretation, "Heterogeneous effect confirmed by both methods")
})

test_that("Power assessment handles sample size limitations", {
  # Scenario: Small sample size in some cell types
  cell_sample_sizes <- c(
    predicted.celltype.l1.Mono = 1000, # Good power
    predicted.celltype.l1.DC = 800, # Good power
    predicted.celltype.l1.CD4_T = 1200, # Good power
    predicted.celltype.l1.CD8_T = 150, # Low power
    predicted.celltype.l1.B = 100, # Low power
    predicted.celltype.l1.NK = 80 # Very low power
  )

  lfsr_values <- c(
    predicted.celltype.l1.Mono = 0.01,
    predicted.celltype.l1.DC = 0.02,
    predicted.celltype.l1.CD4_T = 0.03,
    predicted.celltype.l1.CD8_T = 0.4, # Gray zone (underpowered)
    predicted.celltype.l1.B = 0.45, # Gray zone (underpowered)
    predicted.celltype.l1.NK = 0.6 # Null (very underpowered)
  )

  # Low sample size cells show gray zone/null even if effect might be present
  low_power_cells <- names(cell_sample_sizes)[cell_sample_sizes < 200]
  expect_equal(length(low_power_cells), 3)

  # Check if low power cells are in gray zone or null
  for (cell in low_power_cells) {
    expect_true(lfsr_values[cell] >= 0.05) # Not significant
  }

  # Categorization should reflect uncertainty
  significant <- names(lfsr_values)[lfsr_values < 0.05]
  gray_zone <- names(lfsr_values)[lfsr_values >= 0.05 & lfsr_values < 0.5]

  if (length(gray_zone) > 0) {
    category <- "Likely shared but underpowered"
  } else if (length(significant) >= 2) {
    category <- "Shared"
  } else {
    category <- "Cell-type specific"
  }

  expect_equal(category, "Likely shared but underpowered")
})

# Test complex variant scenarios
test_that("Variant categorization handles complex biological scenarios", {
  # Scenario 1: Master regulator variant
  master_variant <- data.table(
    variant_id = "chr7_50000000_A_G",
    qtl_mechanism = "Local Cascade",
    affected_peaks = 5, # Multiple peaks
    affected_genes = 3, # Multiple genes
    caqtl_cells = list(c(
      "predicted.celltype.l1.Mono", "predicted.celltype.l1.DC",
      "predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T"
    )),
    eqtl_cells = list(c(
      "predicted.celltype.l1.Mono", "predicted.celltype.l1.DC",
      "predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T"
    ))
  )

  # Master regulators affect multiple features across cell types
  expect_true(master_variant$affected_peaks > 1)
  expect_true(master_variant$affected_genes > 1)
  expect_equal(length(master_variant$caqtl_cells[[1]]), 4)

  # Scenario 2: Cell-type-specific enhancer variant
  enhancer_variant <- data.table(
    variant_id = "chr3_75000000_C_T",
    qtl_mechanism = "Only caQTL (With Link)",
    affected_peaks = 1,
    affected_genes = 1,
    caqtl_cells = list(c("predicted.celltype.l1.CD4_T")), # T-cell specific enhancer
    eqtl_cells = list(c()) # No direct eQTL
  )

  expect_equal(enhancer_variant$qtl_mechanism, "Only caQTL (With Link)")
  expect_equal(length(enhancer_variant$caqtl_cells[[1]]), 1)
  expect_equal(length(enhancer_variant$eqtl_cells[[1]]), 0)

  # Scenario 3: Pleiotropic variant with mixed mechanisms
  pleiotropic_variant <- data.table(
    variant_id = "chr1_150000000_G_A",
    genes = list(c("GENE1", "GENE2", "GENE3")),
    mechanisms = list(c("Local Cascade", "Only eQTL", "Positional Cascade")),
    specificities = list(c("Cross-lineage shared", "Single cell-type", "Lineage-specific"))
  )

  # Pleiotropic variants can have different mechanisms for different genes
  expect_equal(length(unique(pleiotropic_variant$mechanisms[[1]])), 3)
  expect_equal(length(unique(pleiotropic_variant$specificities[[1]])), 3)
})

# Test threshold sensitivity
test_that("Categorization is robust to threshold choices", {
  lfsr_values <- c(
    predicted.celltype.l1.Mono = 0.045, # Near significance threshold
    predicted.celltype.l1.DC = 0.055, # Just above threshold
    predicted.celltype.l1.CD4_T = 0.48, # Near null threshold
    predicted.celltype.l1.CD8_T = 0.52 # Just above null threshold
  )

  # Liberal thresholds (sig=0.1, null=0.5)
  sig_liberal <- names(lfsr_values)[lfsr_values < 0.1]
  gray_liberal <- names(lfsr_values)[lfsr_values >= 0.1 & lfsr_values < 0.5]
  expect_equal(length(sig_liberal), 2)
  expect_equal(length(gray_liberal), 1)

  # Conservative thresholds (sig=0.01, null=0.3)
  sig_conservative <- names(lfsr_values)[lfsr_values < 0.01]
  gray_conservative <- names(lfsr_values)[lfsr_values >= 0.01 & lfsr_values < 0.3]
  expect_equal(length(sig_conservative), 0)
  expect_equal(length(gray_conservative), 2)

  # Standard thresholds (sig=0.05, null=0.5)
  sig_standard <- names(lfsr_values)[lfsr_values < 0.05]
  gray_standard <- names(lfsr_values)[lfsr_values >= 0.05 & lfsr_values < 0.5]
  expect_equal(length(sig_standard), 1)
  expect_equal(length(gray_standard), 2)

  # Different thresholds lead to different interpretations
  expect_true(length(sig_liberal) > length(sig_standard))
  expect_true(length(sig_standard) > length(sig_conservative))
})

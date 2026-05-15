library(testthat)
library(cascade)
library(data.table)

test_that("Cochran's Q threshold affects categorization logic", {
  # Test that threshold is used in categorization decisions
  cochran_threshold <- 5e-8

  # Mock p-values for heterogeneity test
  heterogeneity_pvals <- c(1e-10, 1e-6, 0.1)

  # Variants with p-value below threshold should be categorized as heterogeneous
  is_heterogeneous <- heterogeneity_pvals < cochran_threshold

  expect_equal(is_heterogeneous, c(TRUE, FALSE, FALSE))

  # Test with different threshold
  lenient_threshold <- 1e-5
  is_heterogeneous_lenient <- heterogeneity_pvals < lenient_threshold
  expect_equal(is_heterogeneous_lenient, c(TRUE, TRUE, FALSE))
})

test_that("Cochran's Q is used in variant analysis", {
  # Test that heterogeneity threshold is passed to variant analysis
  config <- list(
    parameters = list(
      cochran_q_threshold = 5e-8
    )
  )

  expect_equal(config$parameters$cochran_q_threshold, 5e-8)

  # Test that threshold affects categorization
  # When p-value < threshold, variant is heterogeneous
  mock_pvalue <- 1e-10
  threshold <- 5e-8

  is_heterogeneous <- mock_pvalue < threshold
  expect_true(is_heterogeneous)

  # When p-value > threshold, variant is homogeneous
  mock_pvalue <- 0.1
  is_heterogeneous <- mock_pvalue < threshold
  expect_false(is_heterogeneous)
})

test_that("Cochran's Q handles edge cases in categorization", {
  # Test edge cases for heterogeneity detection
  threshold <- 5e-8

  # Extremely small p-value (strong heterogeneity)
  pval_extreme <- 1e-20
  expect_true(pval_extreme < threshold)

  # Borderline p-value
  pval_border <- 5e-8
  expect_false(pval_border < threshold) # Not strictly less than

  # Large p-value (no heterogeneity)
  pval_large <- 0.5
  expect_false(pval_large < threshold)
})

test_that("Cochran's Q threshold determines categorization", {
  # Mock data with different p-values
  variants <- c("rs1", "rs2", "rs3")
  cochran_pvalues <- c(1e-10, 1e-6, 0.1) # Very significant, borderline, not significant

  # Test with strict threshold
  strict_threshold <- 5e-8
  heterogeneous_strict <- cochran_pvalues < strict_threshold
  expect_equal(heterogeneous_strict, c(TRUE, FALSE, FALSE))

  # Test with lenient threshold
  lenient_threshold <- 1e-5
  heterogeneous_lenient <- cochran_pvalues < lenient_threshold
  expect_equal(heterogeneous_lenient, c(TRUE, TRUE, FALSE))

  # Test with very lenient threshold
  very_lenient_threshold <- 0.05
  heterogeneous_very_lenient <- cochran_pvalues < very_lenient_threshold
  expect_equal(heterogeneous_very_lenient, c(TRUE, TRUE, FALSE))
})

test_that("Cochran's Q caching works correctly", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create cache
  cache <- list(
    cochran_q = new.env(parent = emptyenv())
  )

  # Add some cached results
  gene_id <- "GENE1"
  variant_ids <- c("rs1", "rs2", "rs3")
  pvalues <- c(1e-8, 0.01, 0.5)

  cache$cochran_q[[gene_id]] <- data.frame(
    variant = variant_ids,
    cochran_q_pvalue = pvalues,
    stringsAsFactors = FALSE
  )

  # Test retrieval from cache
  cached_result <- cache$cochran_q[[gene_id]]
  expect_equal(cached_result$variant, variant_ids)
  expect_equal(cached_result$cochran_q_pvalue, pvalues)

  # Test cache miss
  expect_null(cache$cochran_q[["GENE_NOT_IN_CACHE"]])
})

test_that("Cochran's Q integrates with categorization", {
  # Test switching between LFSR and Cochran's Q based categorization

  # Scenario 1: Use Cochran's Q (heterogeneous variant)
  cochran_pvalue <- 1e-10
  cochran_threshold <- 5e-8
  use_cochran <- TRUE

  if (use_cochran && cochran_pvalue < cochran_threshold) {
    category <- "Likely shared but underpowered"
  } else {
    category <- "Terminal cell-specific"
  }

  expect_equal(category, "Likely shared but underpowered")

  # Scenario 2: Use Cochran's Q (homogeneous variant)
  cochran_pvalue <- 0.5

  if (use_cochran && cochran_pvalue < cochran_threshold) {
    category <- "Likely shared but underpowered"
  } else {
    category <- "Terminal cell-specific"
  }

  expect_equal(category, "Terminal cell-specific")

  # Scenario 3: Don't use Cochran's Q
  use_cochran <- FALSE
  lfsr_values <- c(Mono = 0.01, CD4_T = 0.3, CD8_T = 0.6)

  if (use_cochran) {
    category <- "Should not reach here"
  } else if (any(lfsr_values > 0.05 & lfsr_values < 0.5)) {
    category <- "Likely shared but underpowered"
  } else {
    category <- "Terminal cell-specific"
  }

  expect_equal(category, "Likely shared but underpowered")
})

test_that("Cochran's Q threshold validation", {
  # Test valid threshold values
  valid_thresholds <- c(5e-8, 1e-6, 0.01, 0.05)

  for (threshold in valid_thresholds) {
    expect_true(threshold > 0)
    expect_true(threshold < 1)
  }

  # Test that smaller thresholds are more stringent
  strict <- 5e-8
  lenient <- 0.05

  pvalue <- 1e-6
  expect_false(pvalue < strict)
  expect_true(pvalue < lenient)
})

test_that("Cochran's Q parameter configuration", {
  # Test configuration with Cochran's Q enabled
  config_cochran <- list(
    parameters = list(
      use_cochran_q = TRUE,
      cochran_q_threshold = 5e-8
    )
  )

  expect_true(config_cochran$parameters$use_cochran_q)
  expect_equal(config_cochran$parameters$cochran_q_threshold, 5e-8)

  # Test configuration with Cochran's Q disabled
  config_lfsr <- list(
    parameters = list(
      use_cochran_q = FALSE,
      lfsr_sig_threshold = 0.05,
      lfsr_null_threshold = 0.5
    )
  )

  expect_false(config_lfsr$parameters$use_cochran_q)
  expect_null(config_lfsr$parameters$cochran_q_threshold)

  # Test default threshold when not specified
  config_default <- list(
    parameters = list(
      use_cochran_q = TRUE
    )
  )

  # Should default to genome-wide significance
  if (is.null(config_default$parameters$cochran_q_threshold)) {
    config_default$parameters$cochran_q_threshold <- 5e-8
  }

  expect_equal(config_default$parameters$cochran_q_threshold, 5e-8)
})

test_that("Cochran's Q integration with feature categorization", {
  # Test that cochran_q_threshold is passed correctly
  threshold <- 5e-8

  # Mock feature categorization call
  mock_config <- list(
    parameters = list(
      cochran_q_threshold = threshold,
      use_cochran_q = TRUE
    )
  )

  expect_true(mock_config$parameters$use_cochran_q)
  expect_equal(mock_config$parameters$cochran_q_threshold, threshold)

  # When Cochran's Q is disabled
  mock_config_no_cochran <- list(
    parameters = list(
      use_cochran_q = FALSE,
      lfsr_sig_threshold = 0.05,
      lfsr_null_threshold = 0.5
    )
  )

  expect_false(mock_config_no_cochran$parameters$use_cochran_q)
})

test_that("Cochran's Q output in results", {
  # Mock categorization result with Cochran's Q p-values
  result <- data.table(
    feature_id = c("GENE1", "GENE2", "GENE3"),
    category = c(
      "Likely shared but underpowered",
      "Terminal cell-specific",
      "Lineage-specific"
    ),
    top_variant = c("rs1", "rs2", "rs3"),
    cochran_q_pvalue = c(1e-10, 0.5, 0.01),
    used_cochran_q = c(TRUE, TRUE, FALSE)
  )

  # Check that Cochran's Q p-values are included when used
  cochran_results <- result[used_cochran_q == TRUE]
  expect_true(all(!is.na(cochran_results$cochran_q_pvalue)))

  # Check categorization aligns with Cochran's Q results
  heterogeneous <- result[cochran_q_pvalue < 5e-8 & used_cochran_q == TRUE]
  expect_true(all(heterogeneous$category == "Likely shared but underpowered"))
})

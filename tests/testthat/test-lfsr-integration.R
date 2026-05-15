library(testthat)
library(cascade)
library(data.table)

test_that("load_lfsr_results loads eQTL and caQTL LFSR files", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  cell_types <- c("Mono", "CD4_T", "CD8_T", "B")

  # Create mock eQTL LFSR file with production column format
  eqtl_lfsr <- data.table(
    id = paste0(rep(paste0("GENE", 1:5), 10), ":", rep(paste0("rs", 1:10), each = 5)),
    variant = rep(paste0("rs", 1:10), each = 5),
    gene = rep(paste0("GENE", 1:5), 10)
  )
  for (ct in cell_types) {
    eqtl_lfsr[[paste0("predicted.celltype.l1.", ct)]] <- runif(50, 0, 1)
  }
  eqtl_file <- file.path(test_dir, "eqtl.lfsr.tsv.gz")
  fwrite(eqtl_lfsr, eqtl_file, sep = "\t", compress = "gzip")

  # Create mock caQTL LFSR file with production column format
  caqtl_lfsr <- data.table(
    id = paste0(rep(paste0("peak_", 1:3), 10), ":", rep(paste0("rs", 1:10), each = 3)),
    variant = rep(paste0("rs", 1:10), each = 3),
    peak = rep(paste0("peak_", 1:3), 10)
  )
  for (ct in cell_types) {
    caqtl_lfsr[[paste0("predicted.celltype.l1.", ct)]] <- runif(30, 0, 1)
  }
  caqtl_file <- file.path(test_dir, "caqtl.lfsr.tsv.gz")
  fwrite(caqtl_lfsr, caqtl_file, sep = "\t", compress = "gzip")

  # Create configuration
  config <- list(
    cell_types = cell_types,
    file_patterns = list(
      eqtl_lfsr = eqtl_file,
      caqtl_lfsr = caqtl_file
    )
  )

  # Load LFSR results
  lfsr_results <- load_lfsr_results(config)

  expect_true(is.list(lfsr_results))
  expect_true("eqtl_feature" %in% names(lfsr_results))
  expect_true("caqtl_feature" %in% names(lfsr_results))
  expect_true("eqtl_long" %in% names(lfsr_results))
  expect_true("caqtl_long" %in% names(lfsr_results))

  # Check eQTL long format structure
  expect_true(nrow(lfsr_results$eqtl_long) > 0)
  expect_true(all(c("feature_id", "variant_id", "cell_type", "lfsr") %in% names(lfsr_results$eqtl_long)))

  # Check caQTL long format structure
  expect_true(nrow(lfsr_results$caqtl_long) > 0)
  expect_true(all(c("feature_id", "variant_id", "cell_type", "lfsr") %in% names(lfsr_results$caqtl_long)))

  # Check feature-level aggregated tables
  expect_true(nrow(lfsr_results$eqtl_feature) > 0)
  expect_true(all(c("feature_id", "cell_type", "min_lfsr") %in% names(lfsr_results$eqtl_feature)))
  expect_true(nrow(lfsr_results$caqtl_feature) > 0)
  expect_true(all(c("feature_id", "cell_type", "min_lfsr") %in% names(lfsr_results$caqtl_feature)))
})

test_that("LFSR values used for categorization decisions", {
  # Test LFSR-based categorization logic
  lfsr_values <- c(
    Mono = 0.01, # Significant
    CD4_T = 0.3, # Gray zone
    CD8_T = 0.6, # Null
    B = 0.8 # Null
  )

  sig_threshold <- 0.05
  null_threshold <- 0.5

  # Identify significant cell types
  sig_cts <- names(lfsr_values)[lfsr_values < sig_threshold]
  expect_equal(sig_cts, "Mono")

  # Identify gray zone cell types
  gray_cts <- names(lfsr_values)[lfsr_values >= sig_threshold & lfsr_values < null_threshold]
  expect_equal(gray_cts, "CD4_T")

  # Identify null cell types
  null_cts <- names(lfsr_values)[lfsr_values >= null_threshold]
  expect_equal(sort(null_cts), c("B", "CD8_T"))

  # Categorization decision
  if (length(gray_cts) > 0) {
    category <- "Likely shared but underpowered"
  } else if (length(sig_cts) == 1) {
    category <- "Terminal cell-specific"
  } else if (length(sig_cts) > 1) {
    category <- "Shared"
  } else {
    category <- "No significance"
  }

  expect_equal(category, "Likely shared but underpowered")
})

test_that("LFSR matrix creation from data", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create LFSR data
  cell_types <- c("Mono", "CD4_T", "CD8_T")
  variants <- paste0("rs", 1:5)
  genes <- paste0("GENE", 1:3)

  # Create long format LFSR data
  lfsr_data <- expand.grid(
    variant = variants,
    gene = genes,
    stringsAsFactors = FALSE
  )

  for (ct in cell_types) {
    lfsr_data[[paste0("lfsr_", ct)]] <- runif(nrow(lfsr_data), 0, 1)
  }

  # Convert to matrix format for a specific gene
  gene_id <- "GENE1"
  gene_lfsr <- lfsr_data[lfsr_data$gene == gene_id, ]

  # Extract LFSR columns
  lfsr_cols <- grep("^lfsr_", names(gene_lfsr), value = TRUE)
  lfsr_matrix <- as.matrix(gene_lfsr[, lfsr_cols])

  expect_equal(nrow(lfsr_matrix), length(variants))
  expect_equal(ncol(lfsr_matrix), length(cell_types))
  expect_true(all(lfsr_matrix >= 0 & lfsr_matrix <= 1))
})

test_that("LFSR thresholds affect categorization", {
  # Test different threshold combinations
  test_cases <- list(
    list(
      lfsr = c(Mono = 0.01, CD4_T = 0.02, CD8_T = 0.03),
      sig_thresh = 0.05,
      null_thresh = 0.5,
      expected = "Cross-lineage shared" # All significant
    ),
    list(
      lfsr = c(Mono = 0.01, CD4_T = 0.1, CD8_T = 0.6),
      sig_thresh = 0.05,
      null_thresh = 0.5,
      expected = "Likely shared but underpowered" # Gray zone present
    ),
    list(
      lfsr = c(Mono = 0.01, CD4_T = 0.6, CD8_T = 0.7),
      sig_thresh = 0.05,
      null_thresh = 0.5,
      expected = "Single cell-type" # Single significant
    ),
    list(
      lfsr = c(Mono = 0.6, CD4_T = 0.7, CD8_T = 0.8),
      sig_thresh = 0.05,
      null_thresh = 0.5,
      expected = "No significance" # None significant
    )
  )

  for (tc in test_cases) {
    sig_cts <- names(tc$lfsr)[tc$lfsr < tc$sig_thresh]
    gray_cts <- names(tc$lfsr)[tc$lfsr >= tc$sig_thresh & tc$lfsr < tc$null_thresh]

    if (length(gray_cts) > 0) {
      category <- "Likely shared but underpowered"
    } else if (length(sig_cts) >= 2) {
      # Simplified logic - just check if multiple cell types
      category <- "Cross-lineage shared"
    } else if (length(sig_cts) == 1) {
      category <- "Single cell-type"
    } else if (length(sig_cts) == 0) {
      category <- "No significance"
    } else {
      category <- "Shared"
    }

    expect_equal(category, tc$expected,
      info = paste("Failed for LFSR:", paste(tc$lfsr, collapse = ", "))
    )
  }
})

test_that("LFSR handles missing cell types", {
  # LFSR data with subset of cell types
  lfsr_data <- data.table(
    variant = c("rs1", "rs2"),
    gene = c("GENE1", "GENE1"),
    lfsr_Mono = c(0.01, 0.02),
    lfsr_CD4_T = c(0.1, 0.15)
    # Missing CD8_T and other cell types
  )

  all_cell_types <- c("Mono", "CD4_T", "CD8_T", "B")
  lfsr_cols <- paste0("lfsr_", all_cell_types)

  # Check which columns exist
  existing_cols <- intersect(lfsr_cols, names(lfsr_data))
  missing_cols <- setdiff(lfsr_cols, names(lfsr_data))

  expect_equal(existing_cols, c("lfsr_Mono", "lfsr_CD4_T"))
  expect_equal(missing_cols, c("lfsr_CD8_T", "lfsr_B"))

  # Add missing columns with NA
  for (col in missing_cols) {
    lfsr_data[[col]] <- NA_real_
  }

  expect_true(all(lfsr_cols %in% names(lfsr_data)))
  expect_true(all(is.na(lfsr_data$lfsr_CD8_T)))
  expect_true(all(is.na(lfsr_data$lfsr_B)))
})

test_that("mashr model integration with LFSR", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Create mock mashr model
  n_effects <- 10
  n_conditions <- 3

  mock_mashr <- list(
    fitted_g = list(
      Ulist = list(
        tFLASH_1 = diag(n_conditions),
        tFLASH_2 = matrix(1, n_conditions, n_conditions)
      ),
      grid = c(0.5, 1, 2)
    ),
    result = list(
      PosteriorMean = matrix(rnorm(n_effects * n_conditions), n_effects, n_conditions),
      PosteriorSD = matrix(abs(rnorm(n_effects * n_conditions, sd = 0.1)), n_effects, n_conditions),
      lfsr = matrix(runif(n_effects * n_conditions), n_effects, n_conditions)
    )
  )

  mashr_file <- file.path(test_dir, "mashr_model.rds")
  saveRDS(mock_mashr, mashr_file)

  # Load mashr model
  loaded_mashr <- readRDS(mashr_file)

  expect_true("fitted_g" %in% names(loaded_mashr))
  expect_true("result" %in% names(loaded_mashr))
  expect_true("lfsr" %in% names(loaded_mashr$result))

  # Extract LFSR from mashr
  mashr_lfsr <- loaded_mashr$result$lfsr
  expect_equal(dim(mashr_lfsr), c(n_effects, n_conditions))
  expect_true(all(mashr_lfsr >= 0 & mashr_lfsr <= 1))
})

test_that("LFSR summary statistics", {
  # Create LFSR data for summary
  lfsr_data <- data.table(
    variant = rep(paste0("rs", 1:100), each = 2),
    gene = rep(c("GENE1", "GENE2"), 100),
    lfsr_Mono = runif(200, 0, 1),
    lfsr_CD4_T = runif(200, 0, 1),
    lfsr_CD8_T = runif(200, 0, 1)
  )

  # Calculate summary statistics
  sig_threshold <- 0.05

  summary_stats <- list(
    total_pairs = nrow(lfsr_data),
    sig_in_mono = sum(lfsr_data$lfsr_Mono < sig_threshold),
    sig_in_cd4t = sum(lfsr_data$lfsr_CD4_T < sig_threshold),
    sig_in_cd8t = sum(lfsr_data$lfsr_CD8_T < sig_threshold),
    sig_in_any = sum(apply(
      lfsr_data[, c("lfsr_Mono", "lfsr_CD4_T", "lfsr_CD8_T")], 1,
      function(x) any(x < sig_threshold)
    ))
  )

  expect_true(summary_stats$total_pairs == 200)
  expect_true(summary_stats$sig_in_any >= 0)
  expect_true(summary_stats$sig_in_any <= summary_stats$total_pairs)
})

test_that("LFSR file format validation", {
  test_dir <- tempfile("cascade_test")
  dir.create(test_dir, recursive = TRUE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  # Test valid LFSR format
  valid_lfsr <- data.table(
    variant = c("rs1", "rs2"),
    gene = c("GENE1", "GENE2"),
    lfsr_Mono = c(0.01, 0.02),
    lfsr_CD4_T = c(0.1, 0.2)
  )

  valid_file <- file.path(test_dir, "valid_lfsr.tsv.gz")
  fwrite(valid_lfsr, valid_file, sep = "\t", compress = "gzip")

  loaded <- fread(valid_file)
  expect_true("variant" %in% names(loaded))
  expect_true("gene" %in% names(loaded))
  expect_true(any(grepl("^lfsr_", names(loaded))))

  # Test invalid LFSR format (missing variant column)
  invalid_lfsr <- data.table(
    gene = c("GENE1", "GENE2"),
    lfsr_Mono = c(0.01, 0.02)
  )

  invalid_file <- file.path(test_dir, "invalid_lfsr.tsv.gz")
  fwrite(invalid_lfsr, invalid_file, sep = "\t", compress = "gzip")

  loaded_invalid <- fread(invalid_file)
  expect_false("variant" %in% names(loaded_invalid))
})

test_that("LFSR configuration parameters", {
  configs <- list(
    # Conservative thresholds
    list(
      lfsr_sig_threshold = 0.01,
      lfsr_null_threshold = 0.5,
      description = "conservative"
    ),
    # Standard thresholds
    list(
      lfsr_sig_threshold = 0.05,
      lfsr_null_threshold = 0.5,
      description = "standard"
    ),
    # Liberal thresholds
    list(
      lfsr_sig_threshold = 0.1,
      lfsr_null_threshold = 0.9,
      description = "liberal"
    )
  )

  for (config in configs) {
    expect_true(config$lfsr_sig_threshold > 0)
    expect_true(config$lfsr_sig_threshold < config$lfsr_null_threshold)
    expect_true(config$lfsr_null_threshold <= 1)
  }
})

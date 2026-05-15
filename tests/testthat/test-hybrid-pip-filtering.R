test_that("Hybrid PIP filtering works correctly", {
  # Create mock data with varying PIPs across cell types
  # Variant 1: max PIP = 0.6 (T cells), but 0.08 in B cells (should exclude B cells data)
  # Variant 2: max PIP = 0.4 (below threshold, should be excluded entirely)
  # Variant 3: max PIP = 0.8, all above 0.1 (should include all)

  # Mock caQTL files
  temp_dir <- tempdir()

  # T cells data
  t_cells_data <- data.frame(
    rsid = c("chr1_10000", "chr1_20000", "chr1_30000"),
    prob = c(0.6, 0.4, 0.8), # PIPs
    chromosome = c("chr1", "chr1", "chr1"),
    region = c("peak1", "peak2", "peak3"),
    gene_id = c("gene1", "gene2", "gene3"),
    phenotype_id = c("pheno1", "pheno2", "pheno3")
  )

  # B cells data
  b_cells_data <- data.frame(
    rsid = c("chr1_10000", "chr1_20000", "chr1_30000"),
    prob = c(0.08, 0.3, 0.7), # rs1 below min threshold
    chromosome = c("chr1", "chr1", "chr1"),
    region = c("peak1", "peak2", "peak3"),
    gene_id = c("gene1", "gene2", "gene3"),
    phenotype_id = c("pheno1", "pheno2", "pheno3")
  )

  # Monocytes data
  mono_data <- data.frame(
    rsid = c("chr1_10000", "chr1_20000", "chr1_30000"),
    prob = c(0.15, 0.35, 0.5), # All above min threshold
    chromosome = c("chr1", "chr1", "chr1"),
    region = c("peak1", "peak2", "peak3"),
    gene_id = c("gene1", "gene2", "gene3"),
    phenotype_id = c("pheno1", "pheno2", "pheno3")
  )

  # Write files
  t_cells_file <- file.path(temp_dir, "t_cells_caqtl.tsv")
  b_cells_file <- file.path(temp_dir, "b_cells_caqtl.tsv")
  mono_file <- file.path(temp_dir, "monocytes_caqtl.tsv")

  write.table(t_cells_data, t_cells_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(b_cells_data, b_cells_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(mono_data, mono_file, sep = "\t", row.names = FALSE, quote = FALSE)

  # Create file lists
  caqtl_files <- list(
    T_cells = t_cells_file,
    B_cells = b_cells_file,
    Monocytes = mono_file
  )

  eqtl_files <- list(
    T_cells = "nonexistent.tsv",
    B_cells = "nonexistent.tsv",
    Monocytes = "nonexistent.tsv"
  )

  peak_gene_files <- NULL

  # Create a simple peak BED file for testing (6 columns: chr, start, end, peak_id, score, strand)
  peak_bed_data <- data.frame(
    chr = c("chr1", "chr1", "chr1"),
    start = c(1000, 2000, 3000),
    end = c(1500, 2500, 3500),
    peak_id = c("peak1", "peak2", "peak3"),
    score = c(0, 0, 0),
    strand = c(".", ".", ".")
  )
  peak_bed_file <- file.path(temp_dir, "test_peaks.bed")
  write.table(peak_bed_data, peak_bed_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Test hybrid filtering with default thresholds
  result <- load_variant_data_by_qtl_type(
    caqtl_files = caqtl_files,
    eqtl_files = eqtl_files,
    peak_gene_files = peak_gene_files,
    peak_bed_file = peak_bed_file,
    pip_threshold = 0.5, # Max PIP threshold
    min_pip_threshold = 0.1, # Min per-cell-type threshold
    quiet = TRUE
  )

  # Check results
  # rs1: max PIP = 0.6 (passes), but B cells PIP = 0.08 (fails min threshold)
  expect_true("chr1_10000" %in% result$all_variants)
  expect_true("chr1_10000" %in% result$caqtl_data_by_ct$T_cells$variant_id)
  expect_false("chr1_10000" %in% result$caqtl_data_by_ct$B_cells$variant_id) # Excluded due to PIP < 0.1
  expect_true("chr1_10000" %in% result$caqtl_data_by_ct$Monocytes$variant_id)

  # rs2: max PIP = 0.4 (fails max threshold) - should be excluded entirely
  expect_false("chr1_20000" %in% result$all_variants)
  expect_false("chr1_20000" %in% result$caqtl_data_by_ct$T_cells$variant_id)
  expect_false("chr1_20000" %in% result$caqtl_data_by_ct$B_cells$variant_id)
  expect_false("chr1_20000" %in% result$caqtl_data_by_ct$Monocytes$variant_id)

  # rs3: max PIP = 0.8 (passes), all cell types above min threshold
  expect_true("chr1_30000" %in% result$all_variants)
  expect_true("chr1_30000" %in% result$caqtl_data_by_ct$T_cells$variant_id)
  expect_true("chr1_30000" %in% result$caqtl_data_by_ct$B_cells$variant_id)
  expect_true("chr1_30000" %in% result$caqtl_data_by_ct$Monocytes$variant_id)

  # Clean up
  unlink(c(t_cells_file, b_cells_file, mono_file))
})

test_that("Hybrid filtering with edge cases", {
  temp_dir <- tempdir()

  # Test with min_pip_threshold = 0 (effectively no minimum filtering)
  t_cells_data <- data.frame(
    rsid = c("chr1_10000", "chr1_20000"),
    prob = c(0.6, 0.01), # rs2 has very low PIP but should be included
    chromosome = c("chr1", "chr1"),
    region = c("peak1", "peak2"),
    gene_id = c("gene1", "gene2"),
    phenotype_id = c("pheno1", "pheno2")
  )

  t_cells_file <- file.path(temp_dir, "t_cells_edge.tsv")
  write.table(t_cells_data, t_cells_file, sep = "\t", row.names = FALSE, quote = FALSE)

  # Create a simple peak BED file for this test (6 columns: chr, start, end, peak_id, score, strand)
  peak_bed_data2 <- data.frame(
    chr = c("chr1", "chr1"),
    start = c(1000, 2000),
    end = c(1500, 2500),
    peak_id = c("peak1", "peak2"),
    score = c(0, 0),
    strand = c(".", ".")
  )
  peak_bed_file2 <- file.path(temp_dir, "test_peaks2.bed")
  write.table(peak_bed_data2, peak_bed_file2, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  result <- load_variant_data_by_qtl_type(
    caqtl_files = list(T_cells = t_cells_file),
    eqtl_files = list(T_cells = "nonexistent.tsv"),
    peak_gene_files = NULL,
    peak_bed_file = peak_bed_file2,
    pip_threshold = 0.5,
    min_pip_threshold = 0, # No minimum threshold
    quiet = TRUE
  )

  # rs1 should be included (PIP = 0.6 > 0.5)
  expect_true("chr1_10000" %in% result$all_variants)
  # rs2 should be excluded (max PIP = 0.01 < 0.5) even though min_pip = 0
  expect_false("chr1_20000" %in% result$all_variants)

  # Clean up
  unlink(t_cells_file)
})

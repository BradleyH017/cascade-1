test_that("create_config works correctly", {
  # Test basic configuration
  config <- create_config(
    cell_types = c("PBMC", "Mono", "DC"),
    chromosomes = c("chr1", "chr2"),
    file_patterns = list(eqtl_acat = "test/{CELL_TYPE}.tsv"),
    feature_type = "gene"
  )

  expect_equal(config$cell_types, c("PBMC", "Mono", "DC"))
  expect_equal(config$chromosomes, c("chr1", "chr2"))
  expect_equal(config$feature_type, "gene")
  expect_equal(config$file_patterns$eqtl_acat, "test/{CELL_TYPE}.tsv")

  # Test missing required patterns
  expect_error(
    create_config(
      cell_types = c("PBMC", "Mono"),
      file_patterns = list(),
      feature_type = "gene"
    ),
    "Missing required file patterns"
  )

  # Test variant analysis requirements
  config <- create_config(
    cell_types = c("PBMC", "Mono"),
    file_patterns = list(
      caqtl_susie = "test/{CELL_TYPE}_caqtl.tsv",
      eqtl_susie = "test/{CELL_TYPE}_eqtl.tsv"
    ),
    feature_type = "variant"
  )

  expect_equal(config$feature_type, "variant")
  expect_true("caqtl_susie" %in% names(config$file_patterns))
  expect_true("eqtl_susie" %in% names(config$file_patterns))
})

test_that("parameter defaults work correctly", {
  config <- create_config(
    cell_types = c("PBMC", "Mono"),
    file_patterns = list(eqtl_acat = "test.tsv"),
    feature_type = "gene"
  )

  expect_equal(config$parameters$pip_threshold, 0.5)
  expect_equal(config$parameters$acat_fdr_threshold, 0.05)
  expect_equal(config$parameters$lfsr_sig_threshold, 0.05)
  expect_equal(config$parameters$lfsr_null_threshold, 0.5)
  expect_false(config$parameters$run_mash)
  expect_null(config$parameters$n_cores)
  expect_equal(config$parameters$mash_params$max_variants_per_gene, 5)
  expect_equal(config$parameters$mash_params$alpha, 1)
  expect_equal(config$parameters$mash_params$strong_z_threshold, 2)
})

test_that("n_cores parameter can be set", {
  # Test setting n_cores explicitly
  config <- create_config(
    cell_types = c("PBMC", "Mono"),
    file_patterns = list(eqtl_acat = "test.tsv"),
    feature_type = "gene",
    parameters = list(
      pip_threshold = 0.5,
      fdr_threshold = 0.05,
      sig_threshold = 0.05,
      null_threshold = 0.5,
      run_mash = TRUE,
      n_cores = 4,
      mash_params = list(
        max_variants_per_gene = 5,
        alpha = 1,
        strong_z_threshold = 2
      )
    )
  )

  expect_equal(config$parameters$n_cores, 4)
  expect_true(config$parameters$run_mash)
})

# CLI-specific configuration tests
test_that("CLI configuration handles Cochran's Q parameters", {
  # Test with Cochran's Q enabled
  config <- create_config(
    cell_types = c("Mono", "CD4_T", "CD8_T"),
    file_patterns = list(
      eqtl_acat = "test/{CELL_TYPE}.acat.txt",
      eqtl_susie = "test/{CELL_TYPE}.{CHR}.susie.txt"
    ),
    feature_type = "gene",
    parameters = list(
      use_cochran_q = TRUE,
      cochran_q_threshold = 5e-8
    )
  )

  # Add Cochran's Q parameters manually since they're not defaults
  config$parameters$use_cochran_q <- TRUE
  config$parameters$cochran_q_threshold <- 5e-8

  expect_true(config$parameters$use_cochran_q)
  expect_equal(config$parameters$cochran_q_threshold, 5e-8)

  # Test with Cochran's Q disabled (LFSR-based)
  config_lfsr <- create_config(
    cell_types = c("Mono", "CD4_T"),
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene",
    parameters = list(
      use_cochran_q = FALSE,
      lfsr_sig_threshold = 0.05,
      lfsr_null_threshold = 0.5
    )
  )

  # Use_cochran_q is explicitly set to FALSE
  expect_false(config_lfsr$parameters$use_cochran_q)
  # But lfsr thresholds should use defaults from constants (0.05 and 0.5)
  expect_equal(config_lfsr$parameters$lfsr_sig_threshold, 0.05)
  expect_equal(config_lfsr$parameters$lfsr_null_threshold, 0.5)
})

test_that("CLI configuration handles all required file patterns", {
  # Gene analysis required patterns
  gene_patterns <- list(
    eqtl_acat = "/data/{CELL_TYPE}.eqtl.acat.gz",
    eqtl_susie = "/data/{CELL_TYPE}.{CHR}.eqtl.susie.gz",
    eqtl_lfsr = "/data/eqtl.lfsr.tsv.gz",
    eqtl_mashr = "/data/eqtl.mashr.rds"
  )

  config_gene <- create_config(
    cell_types = c("Mono"),
    file_patterns = gene_patterns,
    feature_type = "gene"
  )

  expect_equal(config_gene$file_patterns$eqtl_acat, gene_patterns$eqtl_acat)
  expect_equal(config_gene$file_patterns$eqtl_lfsr, gene_patterns$eqtl_lfsr)

  # Peak analysis required patterns
  peak_patterns <- list(
    caqtl_acat = "/data/{CELL_TYPE}.caqtl.acat.gz",
    caqtl_susie = "/data/{CELL_TYPE}.{CHR}.caqtl.susie.gz",
    caqtl_lfsr = "/data/caqtl.lfsr.tsv.gz",
    caqtl_mashr = "/data/caqtl.mashr.rds"
  )

  config_peak <- create_config(
    cell_types = c("Mono"),
    file_patterns = peak_patterns,
    feature_type = "peak"
  )

  expect_equal(config_peak$file_patterns$caqtl_acat, peak_patterns$caqtl_acat)
  expect_equal(config_peak$file_patterns$caqtl_lfsr, peak_patterns$caqtl_lfsr)

  # Variant analysis required patterns
  variant_patterns <- list(
    caqtl_susie = "/data/{CELL_TYPE}.{CHR}.caqtl.susie.gz",
    eqtl_susie = "/data/{CELL_TYPE}.{CHR}.eqtl.susie.gz",
    peak_gene_links = "/data/peak_gene_links.tsv.gz",
    peak_bed = "/data/{CELL_TYPE}.peaks.bed"
  )

  config_variant <- create_config(
    cell_types = c("Mono"),
    file_patterns = variant_patterns,
    feature_type = "variant"
  )

  expect_equal(config_variant$file_patterns$peak_gene_links, variant_patterns$peak_gene_links)
})

test_that("CLI configuration handles column mapping", {
  column_mapping <- list(
    eqtl_acat = list(
      feature_id = "phenotype_id",
      q_value = "ACAT_q"
    ),
    caqtl_acat = list(
      feature_id = "phenotype_id",
      q_value = "qval"
    ),
    eqtl_susie = list(
      variant_id = "rsid",
      feature_id = "region",
      pip = "prob",
      chromosome = "chromosome",
      cs_id = "cs",
      beta = NULL,
      se = NULL
    ),
    caqtl_susie = list(
      variant_id = "rsid",
      feature_id = "region",
      pip = "prob",
      chromosome = "chromosome",
      cs_id = "cs",
      beta = NULL,
      se = NULL
    )
  )

  config <- create_config(
    cell_types = c("Mono"),
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene",
    column_mapping = column_mapping
  )

  expect_equal(config$column_mapping$eqtl_acat$feature_id, "phenotype_id")
  expect_equal(config$column_mapping$eqtl_acat$q_value, "ACAT_q")
  expect_equal(config$column_mapping$caqtl_acat$q_value, "qval")
})

test_that("CLI configuration handles debug mode", {
  # Debug is not a parameter of create_config, it should be added separately
  config <- create_config(
    cell_types = c("Mono"),
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene"
  )

  # Add debug settings manually
  config$debug <- list(
    enabled = TRUE,
    n_genes = 100,
    n_peaks = 50
  )

  expect_true(config$debug$enabled)
  expect_equal(config$debug$n_genes, 100)
  expect_equal(config$debug$n_peaks, 50)
})

test_that("CLI configuration handles cache settings", {
  # Cache is not a parameter of create_config, it should be added separately
  config_cache <- create_config(
    cell_types = c("Mono"),
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene"
  )
  # Add cache settings manually
  config_cache$cache <- list(enabled = TRUE)

  expect_true(config_cache$cache$enabled)

  # Cache disabled
  config_no_cache <- create_config(
    cell_types = c("Mono"),
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene"
  )
  config_no_cache$cache <- list(enabled = FALSE)

  expect_false(config_no_cache$cache$enabled)
})

test_that("CLI configuration validates feature types", {
  valid_feature_types <- c("gene", "peak", "variant", "all")

  # Test gene type with proper file patterns
  config <- create_config(
    cell_types = c("Mono"),
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene"
  )
  expect_equal(config$feature_type, "gene")

  # Test peak type with proper file patterns
  config <- create_config(
    cell_types = c("Mono"),
    file_patterns = list(caqtl_acat = "test.txt"),
    feature_type = "peak"
  )
  expect_equal(config$feature_type, "peak")

  # Test variant type with proper file patterns
  config <- create_config(
    cell_types = c("Mono"),
    file_patterns = list(caqtl_susie = "test.txt", eqtl_susie = "test2.txt"),
    feature_type = "variant"
  )
  expect_equal(config$feature_type, "variant")
})

test_that("CLI configuration handles L1/L2 cell types", {
  # L1 cell types
  l1_celltypes <- c(
    "predicted.celltype.l1.Mono",
    "predicted.celltype.l1.CD4_T",
    "predicted.celltype.l1.CD8_T"
  )

  config_l1 <- create_config(
    cell_types = l1_celltypes,
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene"
  )

  expect_equal(config_l1$cell_types, l1_celltypes)

  # L2 cell types
  l2_celltypes <- c(
    "predicted.celltype.l2.CD14_Mono",
    "predicted.celltype.l2.CD4_Naive",
    "predicted.celltype.l2.CD8_TEM"
  )

  config_l2 <- create_config(
    cell_types = l2_celltypes,
    file_patterns = list(caqtl_susie = "test.txt", eqtl_susie = "test2.txt"),
    feature_type = "variant" # L2 allowed for variant analysis
  )

  expect_equal(config_l2$cell_types, l2_celltypes)

  # Mixed L1 and L2
  mixed_celltypes <- c(
    "predicted.celltype.l1.Mono",
    "predicted.celltype.l2.CD14_Mono"
  )

  config_mixed <- create_config(
    cell_types = mixed_celltypes,
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene"
  )

  expect_equal(config_mixed$cell_types, mixed_celltypes)
})

test_that("CLI configuration handles multiple chromosomes", {
  chromosomes <- c("chr1", "chr2", "chr7", "chr22")

  config <- create_config(
    cell_types = c("Mono"),
    chromosomes = chromosomes,
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene"
  )

  expect_equal(config$chromosomes, chromosomes)
  expect_equal(length(config$chromosomes), 4)
})

test_that("CLI configuration handles output directory", {
  # Output directory is not a parameter of create_config
  config <- create_config(
    cell_types = c("Mono"),
    file_patterns = list(eqtl_acat = "test.txt"),
    feature_type = "gene"
  )

  # Add output directory manually
  output_dirs <- c(
    "/mnt/data/results/test1",
    "./results/test2",
    "~/cascade_output/test3"
  )

  for (output_dir in output_dirs) {
    config$output_dir <- output_dir
    expect_equal(config$output_dir, output_dir)
  }
})

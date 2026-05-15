# test-cache-correctness.R
#
# Tests for CASCADE cache behavior: cached_fread, cache invalidation,
# cache key differentiation, and config hash generation.

test_that("cached_fread with select= returns correct columns", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  data.table::fwrite(data.table::data.table(a = 1:3, b = 4:6, c = 7:9), tmp)

  cache_dir <- file.path(tempdir(), "cache_select_test")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  cache <- memoise::cache_filesystem(cache_dir)
  attr(cache, "cache_dir") <- cache_dir

  result <- cascade:::cached_fread(tmp, cache = cache, select = c("a", "c"))
  expect_equal(names(result), c("a", "c"))
  expect_equal(nrow(result), 3L)
})

test_that("cache hit returns identical data", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  data.table::fwrite(data.table::data.table(x = 1:5, y = letters[1:5]), tmp)


  cache_dir <- file.path(tempdir(), "cache_hit_test")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  cache <- memoise::cache_filesystem(cache_dir)
  attr(cache, "cache_dir") <- cache_dir

  first <- cascade:::cached_fread(tmp, cache = cache)
  second <- cascade:::cached_fread(tmp, cache = cache)
  expect_true(all.equal(first, second))
})

test_that("cache invalidation on mtime change", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  data.table::fwrite(data.table::data.table(val = 1:3), tmp)

  cache_dir <- file.path(tempdir(), "cache_mtime_test")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  cache <- memoise::cache_filesystem(cache_dir)
  attr(cache, "cache_dir") <- cache_dir

  first <- cascade:::cached_fread(tmp, cache = cache)

  # Overwrite with different data and ensure mtime changes
  Sys.sleep(1.1)
  data.table::fwrite(data.table::data.table(val = 10:12), tmp)

  second <- cascade:::cached_fread(tmp, cache = cache)
  expect_equal(second$val, 10:12)
})

test_that("different select= produces different cache entries", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  data.table::fwrite(data.table::data.table(a = 1:3, b = 4:6, c = 7:9), tmp)

  cache_dir <- file.path(tempdir(), "cache_diffselect_test")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  cache <- memoise::cache_filesystem(cache_dir)
  attr(cache, "cache_dir") <- cache_dir

  res_ab <- cascade:::cached_fread(tmp, cache = cache, select = c("a", "b"))
  res_bc <- cascade:::cached_fread(tmp, cache = cache, select = c("b", "c"))

  expect_equal(names(res_ab), c("a", "b"))
  expect_equal(names(res_bc), c("b", "c"))
})

test_that("config hash includes parameters", {
  base_config <- list(
    file_patterns = list(acat = "acat*.tsv"),
    cell_types = c("Mono", "NK"),
    parameters = list(pip_threshold = 0.5),
    debug = FALSE
  )
  hash1 <- cascade:::generate_config_hash(base_config)

  changed_config <- base_config
  changed_config$parameters$pip_threshold <- 0.9
  hash2 <- cascade:::generate_config_hash(changed_config)

  expect_false(hash1 == hash2)
})

# Regression tests for how results print and combine. The bug these guard
# against: output = "tibble" used to keep the report-printing class, so a
# returned summary rendered as a one-pair report - and once several were
# row-bound, printing showed only the first pair while silently hiding the
# rest.

test_that("output = 'tibble' returns something that prints as a tibble", {
  geno <- synthetic_geno()

  cmp <- compare_genotypes(geno, "kid", "mom", output = "tibble")
  expect_false(inherits(cmp, "geno_comparison"))
  expect_s3_class(cmp, "tbl_df")
  expect_output(print(cmp), "n_loci_compared")
  expect_failure(expect_output(print(cmp), "Genotype comparison"))

  trio <- check_trio(geno, "kid", "mom", "dad", output = "tibble")
  expect_false(inherits(trio, "geno_trio"))
  expect_s3_class(trio, "tbl_df")
  expect_failure(expect_output(print(trio), "Trio parentage check"))
})

test_that("locus_detail() still works on tibble-mode results", {
  geno <- synthetic_geno()

  cmp <- compare_genotypes(geno, "kid", "mom", output = "tibble")
  expect_equal(nrow(locus_detail(cmp)), length(get_loci(geno)))

  trio <- check_trio(geno, "kid", "mom", "dad", output = "tibble")
  expect_equal(nrow(locus_detail(trio)), length(get_loci(geno)))
})

test_that("row-binding many comparisons keeps every row visible", {
  geno <- synthetic_geno()
  pairs <- list(c("kid", "mom"), c("kid", "dad"), c("kid", "unrelated1"))
  bound <- do.call(rbind, lapply(pairs, function(p) {
    compare_genotypes(geno, p[1], p[2], output = "tibble")
  }))

  expect_equal(nrow(bound), 3)
  out <- utils::capture.output(print(bound))
  # all three rows must be shown, not just the first pair's report
  # (long IDs get truncated in tibble display, so assert on the row count
  # tibble reports plus a short ID from a later row)
  expect_true(any(grepl("A tibble: 3", out)))
  expect_true(any(grepl("\\bdad\\b", out)))
  expect_false(any(grepl("Genotype comparison", out)))
})

test_that("a multi-row object still prints as a tibble even if classed", {
  # belt-and-braces: even if a geno_comparison-classed object somehow ends up
  # with several rows, the single-pair report must not be used for it
  geno <- synthetic_geno()
  one <- suppressWarnings(utils::capture.output(
    res <- compare_genotypes(geno, "kid", "mom")
  ))
  two <- rbind(res, res)
  expect_s3_class(two, "geno_comparison")
  expect_failure(expect_output(print(two), "Genotype comparison"))
})

test_that("pretty mode still prints the single-pair report", {
  geno <- synthetic_geno()
  expect_output(compare_genotypes(geno, "kid", "mom"), "Genotype comparison")
  expect_output(check_trio(geno, "kid", "mom", "dad"), "Trio parentage check")
})

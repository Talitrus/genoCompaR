test_that("find_duplicate_genotypes flags identical genotypes", {
  geno <- synthetic_geno()
  dupes <- find_duplicate_genotypes(geno, threshold = 0.95, min_loci_compared = 1)
  expect_s3_class(dupes, "geno_duplicates")
  expect_equal(nrow(dupes), 1)
  expect_setequal(c(dupes$sample1, dupes$sample2), c("kid", "dup_a"))
  expect_equal(dupes$pct_match, 100)
})

test_that("lowering the threshold picks up near-duplicates", {
  geno <- synthetic_geno()
  dupes <- find_duplicate_genotypes(geno, threshold = 0.75, min_loci_compared = 1)
  expect_gte(nrow(dupes), 2)
  expect_true(any(round(dupes$pct_match) == 80))
  expect_true(is.unsorted(-dupes$pct_match) == FALSE)
})

test_that("min_loci_compared filters out low-information pairs", {
  geno <- synthetic_geno()
  dupes <- find_duplicate_genotypes(geno, threshold = 0.5, min_loci_compared = 100)
  expect_equal(nrow(dupes), 0)
})

test_that("group columns are carried through when present", {
  geno <- synthetic_geno()
  dupes <- find_duplicate_genotypes(geno, threshold = 0.95, min_loci_compared = 1)
  expect_true(all(c("group1", "group2") %in% names(dupes)))
})

test_that("find_duplicate_genotypes validates inputs", {
  geno <- synthetic_geno()
  expect_error(find_duplicate_genotypes(geno, threshold = 1.5), "threshold")
  expect_error(find_duplicate_genotypes(unclass(geno)), "geno_data")
})

test_that("find_duplicate_genotypes flags well-formed, correctly-thresholded pairs on real data", {
  geno <- read_genotypes(example_data_path())
  dupes <- find_duplicate_genotypes(geno, threshold = 0.95, min_loci_compared = 50)
  expect_gt(nrow(dupes), 0)
  expect_true(all(dupes$sample1 != dupes$sample2))
  expect_true(all(dupes$pct_match >= 95))
  expect_true(all(dupes$n_loci_compared >= 50))
  expect_true(all(dupes$n_match + dupes$n_mismatch == dupes$n_loci_compared))
  expect_false(is.unsorted(rev(dupes$pct_match)))
})

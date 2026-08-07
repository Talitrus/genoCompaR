test_that("compare_genotypes computes match stats and handles a missing locus", {
  geno <- synthetic_geno()
  res <- compare_genotypes(geno, "kid", "mom", output = "tibble")
  expect_equal(res$n_loci_total, 5)
  expect_equal(res$n_loci_compared, 4)
  expect_equal(res$n_loci_unique, 1)
  expect_equal(res$n_match, 1)
  expect_equal(res$n_mismatch, 3)
  expect_equal(res$pct_match, 25)
})

test_that("a true parent-offspring pair has few Mendelian incompatibilities", {
  geno <- synthetic_geno()
  res_mom <- compare_genotypes(geno, "kid", "mom", output = "tibble")
  res_dad <- compare_genotypes(geno, "kid", "dad", output = "tibble")
  expect_equal(res_mom$n_mendelian_incompatible, 1)
  expect_equal(res_dad$n_mendelian_incompatible, 1)
})

test_that("compare_genotypes flags Mendelian incompatibilities bi-directionally", {
  geno <- synthetic_geno()
  res <- compare_genotypes(geno, "kid", "unrelated1", output = "tibble")
  expect_equal(res$n_mendelian_incompatible, 4)

  detail <- locus_detail(res)
  bad <- detail[detail$mendelian_incompatible %in% TRUE, ]
  expect_true(nrow(bad) == 4)
  expect_true(all(!is.na(bad$sample1_alleles_not_in_sample2)))
  expect_true(all(!is.na(bad$sample2_alleles_not_in_sample1)))
})

test_that("mendelian = FALSE skips incompatibility computation", {
  geno <- synthetic_geno()
  res <- compare_genotypes(geno, "kid", "unrelated1", mendelian = FALSE, output = "tibble")
  expect_false("n_mendelian_incompatible" %in% names(res))
  detail <- locus_detail(res)
  expect_true(all(is.na(detail$mendelian_incompatible)))
})

test_that("identical genotypes give 100% match and zero incompatibilities", {
  geno <- synthetic_geno()
  res <- compare_genotypes(geno, "kid", "dup_a", output = "tibble")
  expect_equal(res$pct_match, 100)
  expect_equal(res$n_mendelian_incompatible, 0)
})

test_that("compare_genotypes validates inputs", {
  geno <- synthetic_geno()
  expect_error(compare_genotypes(geno, "kid", "kid"), "different")
  expect_error(compare_genotypes(geno, "kid", "nobody"), "not found")
  expect_error(compare_genotypes(unclass(geno), "kid", "mom"), "geno_data")
})

test_that("pretty output prints a human-readable report and returns invisibly", {
  geno <- synthetic_geno()
  expect_output(res <- compare_genotypes(geno, "kid", "mom"), "Genotype comparison")
  expect_output(compare_genotypes(geno, "kid", "unrelated1"), "Mendelian incompatibilities")
  expect_s3_class(res, "geno_comparison")
})

test_that("compare_genotypes agrees with find_duplicate_genotypes on real data", {
  # Cross-validate against an independently-computed pair rather than
  # asserting specifics about this dataset's biology: pull the top match
  # from find_duplicate_genotypes() and confirm compare_genotypes() (which
  # uses a completely different, non-vectorized code path) reports the same
  # counts for that exact pair.
  geno <- read_genotypes(example_data_path())
  dupes <- find_duplicate_genotypes(geno, threshold = 0.95, min_loci_compared = 50)
  expect_gt(nrow(dupes), 0)

  top <- dupes[1, ]
  res <- compare_genotypes(geno, top$sample1, top$sample2, output = "tibble")
  expect_equal(res$n_loci_compared, top$n_loci_compared)
  expect_equal(res$n_match, top$n_match)
  expect_equal(res$pct_match, top$pct_match)
})

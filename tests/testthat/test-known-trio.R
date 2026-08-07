# Validates check_trio() against a real known-truth parent-offspring trio
# (broodstock reference genotypes vs. a salvage-sample offspring), rather
# than only synthetic data. Deliberately avoids any all-pairs scan across
# either source file (which are ~3,200 and ~1,200 samples respectively) -
# only the three individuals of interest are ever loaded into a geno_data
# object.

test_that("check_trio confirms a real known parent-offspring trio", {
  geno <- .known_trio_geno()

  res <- check_trio(geno, "C260383CVP", "C240716FRH", "C240722FRH", output = "tibble")
  expect_s3_class(res, "tbl_df")
  expect_equal(res$n_loci_compared, 178)
  expect_equal(res$n_compatible, 176)
  expect_equal(res$n_incompatible, 2)
  expect_equal(res$pct_compatible, 176 / 178 * 100)

  # a true trio should be compatible at the overwhelming majority of loci -
  # the couple of misses are consistent with ordinary genotyping error
  expect_gt(res$pct_compatible, 95)
})

test_that("check_trio is direction-agnostic for the real trio (parent order doesn't matter)", {
  geno <- .known_trio_geno()

  a <- check_trio(geno, "C260383CVP", "C240716FRH", "C240722FRH", output = "tibble")
  b <- check_trio(geno, "C260383CVP", "C240722FRH", "C240716FRH", output = "tibble")
  expect_equal(a$n_compatible, b$n_compatible)
  expect_equal(a$n_incompatible, b$n_incompatible)
})

test_that("check_trio discriminates a true parent from a random non-parent", {
  geno <- .known_trio_geno(parent_ids = c("C240716FRH", "C240722FRH", "C240069FRH"))

  true_trio <- check_trio(geno, "C260383CVP", "C240716FRH", "C240722FRH", output = "tibble")
  impostor_trio <- check_trio(geno, "C260383CVP", "C240716FRH", "C240069FRH", output = "tibble")

  expect_gt(true_trio$pct_compatible, impostor_trio$pct_compatible)
  expect_gt(true_trio$pct_compatible - impostor_trio$pct_compatible, 10)
})

test_that("duo-level Mendelian incompatibilities also fail the trio check", {
  geno <- .known_trio_geno()

  duo1 <- compare_genotypes(geno, "C260383CVP", "C240716FRH", output = "tibble")
  duo2 <- compare_genotypes(geno, "C260383CVP", "C240722FRH", output = "tibble")
  trio <- check_trio(geno, "C260383CVP", "C240716FRH", "C240722FRH", output = "tibble")

  # a genuine parent-offspring pair should show very few duo incompatibilities
  expect_lte(duo1$n_mendelian_incompatible, 3)
  expect_lte(duo2$n_mendelian_incompatible, 3)

  # any locus where the offspring shares no allele at all with one parent
  # cannot possibly admit a one-from-each assignment, so it must also be
  # incompatible in the trio result
  duo_bad <- unique(c(
    locus_detail(duo1)$locus[locus_detail(duo1)$mendelian_incompatible %in% TRUE],
    locus_detail(duo2)$locus[locus_detail(duo2)$mendelian_incompatible %in% TRUE]
  ))
  trio_detail <- locus_detail(trio)
  trio_bad <- trio_detail$locus[trio_detail$status == "incompatible"]
  expect_true(all(duo_bad %in% trio_bad))
  expect_gt(length(duo_bad), 0) # guard against a vacuously-true subset check
})

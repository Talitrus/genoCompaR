test_that("locus_missingness counts missing calls per locus", {
  geno <- synthetic_geno()
  lm <- locus_missingness(geno)
  expect_equal(lm$n_missing[lm$locus == "locus4"], 1)
  expect_equal(lm$n_missing[lm$locus == "locus1"], 0)
  expect_equal(lm$pct_missing[lm$locus == "locus4"], 100 / 6)
})

test_that("sample_missingness counts missing loci per sample", {
  geno <- synthetic_geno()
  sm <- sample_missingness(geno)
  expect_equal(sm$n_loci_missing[sm$indiv.ID == "mom"], 1)
  expect_equal(sm$n_loci_missing[sm$indiv.ID == "kid"], 0)
})

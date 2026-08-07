test_that("genotype_match_matrix reproduces compare_genotypes pct_match", {
  geno <- synthetic_geno()
  mm <- genotype_match_matrix(geno, min_loci_compared = 1)
  expect_s3_class(mm, "geno_match_matrix")
  expect_equal(dim(mm), c(6, 6))
  expect_equal(rownames(mm), geno$indiv.ID)

  expect_equal(mm["kid", "dup_a"], 100)
  expect_equal(mm["kid", "dup_b"], 80)
  expect_equal(mm["kid", "unrelated1"], 0)
  expect_true(all(mm == t(mm), na.rm = TRUE))
})

test_that("diagonal is always NA", {
  geno <- synthetic_geno()
  mm <- genotype_match_matrix(geno)
  expect_true(all(is.na(diag(mm))))
})

test_that("min_loci_compared masks low-information cells with NA", {
  geno <- synthetic_geno()
  mm <- genotype_match_matrix(geno, min_loci_compared = 100)
  off_diag <- mm[upper.tri(mm)]
  expect_true(all(is.na(off_diag)))
})

test_that("genotype_match_matrix agrees with compare_genotypes cell-by-cell", {
  geno <- synthetic_geno()
  mm <- genotype_match_matrix(geno, min_loci_compared = 1)
  cmp <- compare_genotypes(geno, "kid", "mom", output = "tibble")
  expect_equal(mm["kid", "mom"], cmp$pct_match)
})

test_that("genotype_match_matrix validates inputs", {
  geno <- synthetic_geno()
  expect_error(genotype_match_matrix(unclass(geno)), "geno_data")
})

test_that("arithmetic drops the match-matrix class and stale attributes", {
  # otherwise `1 - mm / 100` would still print a header describing match
  # percentages, and carry an n_compared attribute that no longer applies
  geno <- synthetic_geno()
  mm <- genotype_match_matrix(geno, min_loci_compared = 1)
  d <- 1 - mm / 100

  expect_false(inherits(d, "geno_match_matrix"))
  expect_null(attr(d, "n_compared"))
  expect_true(is.matrix(d))
  expect_equal(d["kid", "dup_a"], 0) # 100% match -> zero distance
  expect_false(inherits(-mm, "geno_match_matrix")) # unary op too
})

test_that("the documented clustering workflow runs", {
  geno <- read_genotypes(example_data_path())
  mm <- genotype_match_matrix(geno)
  expect_true(any(is.na(mm[upper.tri(mm)]))) # real panels do have NA cells

  d <- stats::as.dist(1 - mm / 100)
  d[is.na(d)] <- 1
  h <- stats::hclust(d)
  expect_s3_class(h, "hclust")
  expect_equal(length(h$order), nrow(geno))
})

test_that("genotype_match_matrix is well-formed on real data", {
  geno <- read_genotypes(example_data_path())
  mm <- genotype_match_matrix(geno, min_loci_compared = 50)
  expect_s3_class(mm, "geno_match_matrix")
  expect_equal(dim(mm), c(nrow(geno), nrow(geno)))
  expect_true(all(is.na(diag(mm))))
  vals <- mm[upper.tri(mm)]
  vals <- vals[!is.na(vals)]
  expect_true(all(vals >= 0 & vals <= 100))
})

test_that("check_trio finds a true parent-offspring trio mostly compatible", {
  geno <- synthetic_geno()
  res <- check_trio(geno, "kid", "mom", "dad", output = "tibble")
  expect_s3_class(res, "tbl_df")
  expect_equal(res$n_loci_compared, 4)
  expect_equal(res$n_loci_unique, 1) # locus4 missing in mom
  expect_equal(res$n_compatible, 3)
  expect_equal(res$n_incompatible, 1) # locus5: kid C/C, neither parent has C
  expect_equal(res$pct_compatible, 75)
})

test_that("check_trio flags a non-parent as mostly incompatible", {
  geno <- synthetic_geno()
  res <- check_trio(geno, "kid", "mom", "unrelated1", output = "tibble")
  expect_equal(res$n_incompatible, 3)
  expect_equal(res$pct_compatible, 25)
})

test_that("check_trio is a stronger test than two independent duo checks", {
  # offspring heterozygous {A,B}; parent1 homozygous {A,A}; parent2 has {A,C}.
  # Both duo checks pass (each shares allele A), but no valid one-allele-
  # from-each assignment exists, so the trio should be incompatible.
  df <- data.frame(
    group = "G1",
    indiv.ID = c("off", "p1", "p2"),
    locusX.1 = c("A", "A", "A"),
    locusX.2 = c("B", "A", "C"),
    stringsAsFactors = FALSE
  )
  geno <- read_genotypes(synthetic_geno_csv(df))

  duo1 <- compare_genotypes(geno, "off", "p1", output = "tibble")
  duo2 <- compare_genotypes(geno, "off", "p2", output = "tibble")
  expect_equal(duo1$n_mendelian_incompatible, 0)
  expect_equal(duo2$n_mendelian_incompatible, 0)

  trio <- check_trio(geno, "off", "p1", "p2", output = "tibble")
  expect_equal(trio$n_incompatible, 1)
  expect_equal(trio$n_compatible, 0)
})

test_that("check_trio treats homozygous offspring correctly", {
  df <- data.frame(
    group = "G1",
    indiv.ID = c("off", "p1", "p2"),
    locusX.1 = c("A", "A", "G"),
    locusX.2 = c("A", "G", "A"),
    stringsAsFactors = FALSE
  )
  geno <- read_genotypes(synthetic_geno_csv(df))
  trio <- check_trio(geno, "off", "p1", "p2", output = "tibble")
  expect_equal(trio$n_compatible, 1)
})

test_that("parent1/parent2 order does not matter", {
  geno <- synthetic_geno()
  a <- check_trio(geno, "kid", "mom", "dad", output = "tibble")
  b <- check_trio(geno, "kid", "dad", "mom", output = "tibble")
  expect_equal(a$n_compatible, b$n_compatible)
  expect_equal(a$n_incompatible, b$n_incompatible)
})

test_that("check_trio validates inputs", {
  geno <- synthetic_geno()
  expect_error(check_trio(geno, "kid", "mom", "kid"), "different")
  expect_error(check_trio(geno, "kid", "mom", "nobody"), "not found")
  expect_error(check_trio(unclass(geno), "kid", "mom", "dad"), "geno_data")
})

test_that("pretty output prints a human-readable report", {
  geno <- synthetic_geno()
  expect_output(res <- check_trio(geno, "kid", "mom", "dad"), "Trio parentage check")
  expect_s3_class(res, "geno_trio")
})

test_that("check_trio is well-formed on real data", {
  geno <- read_genotypes(example_data_path())
  res <- check_trio(geno, "CH-00001_S1", "CH-00002_S2", "CH-00003_S3", output = "tibble")
  expect_s3_class(res, "tbl_df")
  expect_gt(res$n_loci_compared, 0)
  expect_equal(res$n_compatible + res$n_incompatible, res$n_loci_compared)
  expect_equal(res$n_loci_compared + res$n_loci_unique, res$n_loci_total)
})

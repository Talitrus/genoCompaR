test_that("read_genotypes parses a well-formed file", {
  geno <- synthetic_geno()
  expect_s3_class(geno, "geno_data")
  expect_equal(nrow(geno), 6)
  expect_equal(get_loci(geno), paste0("locus", 1:5))
  expect_equal(attr(geno, "id_col"), "indiv.ID")
  expect_equal(attr(geno, "group_col"), "group")
})

test_that("read_genotypes errors on a missing id_col", {
  path <- synthetic_geno_csv()
  expect_error(read_genotypes(path, id_col = "nope"), "id_col")
})

test_that("read_genotypes errors on duplicate sample IDs", {
  df <- synthetic_geno_df()
  df$indiv.ID[2] <- df$indiv.ID[1]
  expect_error(read_genotypes(synthetic_geno_csv(df)), "Duplicate sample IDs")
})

test_that("read_genotypes errors on an unpaired locus column", {
  df <- synthetic_geno_df()
  names(df)[names(df) == "locus5.2"] <- "locus6.2"
  expect_error(read_genotypes(synthetic_geno_csv(df)), "locus5")
})

test_that("read_genotypes errors on an odd number of locus columns", {
  df <- synthetic_geno_df()
  df$stray <- "A"
  expect_error(read_genotypes(synthetic_geno_csv(df)), "even number")
})

test_that("read_genotypes works without a group column", {
  df <- synthetic_geno_df()
  df$group <- NULL
  geno <- read_genotypes(synthetic_geno_csv(df))
  expect_true(is.na(attr(geno, "group_col")))
  expect_equal(nrow(geno), 6)
})

test_that("read_genotypes loads the real example dataset", {
  geno <- read_genotypes(example_data_path())
  expect_s3_class(geno, "geno_data")
  expect_equal(nrow(geno), 288)
  expect_equal(length(get_loci(geno)), 193)
})

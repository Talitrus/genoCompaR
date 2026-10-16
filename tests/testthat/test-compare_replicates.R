test_that("compare_replicates input validation works", {
  geno <- synthetic_geno()

  # Error if not geno_data
  expect_error(compare_replicates(as.data.frame(geno)), "geno_data")

  # Error if < 2 samples
  expect_error(compare_replicates(geno, samples = "kid"), "at least 2 samples")

  # Error if sample IDs not found
  expect_error(compare_replicates(geno, samples = c("kid", "nonexistent")), "Sample\\(s\\) not found")

  # Error if numeric indices out of bounds
  expect_error(compare_replicates(geno, samples = c(1, 99)), "between 1 and")

  # Error if unknown loci requested
  expect_error(compare_replicates(geno, loci = "nonexistent_locus"), "loci were not found")

  # Warning and deduplication if duplicate samples passed
  expect_warning(res <- compare_replicates(geno, samples = c("kid", "dup_a", "kid")), "deduplicated")
  expect_equal(res$samples, c("kid", "dup_a"))
})

test_that("compare_replicates defaults to all samples when samples = NULL", {
  geno <- synthetic_geno()
  sub_geno <- geno[1:3, ]
  res <- compare_replicates(sub_geno, min_loci_compared = 1)

  expect_equal(res$samples, c("mom", "dad", "kid"))
  expect_equal(nrow(res$summary), 3) # 3 choose 2
})

test_that("identical replicates have distance 0 and no discrepancies", {
  geno <- synthetic_geno()
  res <- compare_replicates(geno, samples = c("kid", "dup_a"), min_loci_compared = 1)

  expect_equal(res$summary$n_match, 5)
  expect_equal(res$summary$n_mismatch, 0)
  expect_equal(res$summary$pct_match, 100)
  expect_equal(res$summary$distance, 0)
  expect_equal(nrow(res$discrepancies), 0)
  expect_s3_class(res$distance_matrix, "dist")
  expect_equal(as.matrix(res$distance_matrix)["kid", "dup_a"], 0)
})

test_that("discrepant loci are correctly identified and categorized", {
  df <- data.frame(
    indiv.ID = c("rep1", "rep2", "rep3"),
    # locus1: perfect match across all three
    loc1.1 = c("A", "A", "A"),
    loc1.2 = c("T", "T", "T"),
    # locus2: allelic dropout (rep1 is A/T, rep2 is A/A, rep3 is A/T)
    loc2.1 = c("A", "A", "A"),
    loc2.2 = c("T", "A", "T"),
    # locus3: discordant (rep1 is A/A, rep2 is G/G, rep3 is A/A)
    loc3.1 = c("A", "G", "A"),
    loc3.2 = c("A", "G", "A"),
    # locus4: multiallelic discordant (rep1 is A/T, rep2 is A/C, rep3 is G/T)
    loc4.1 = c("A", "A", "G"),
    loc4.2 = c("T", "C", "T"),
    # locus5: uncalled in one sample, concordant in other two
    loc5.1 = c("C", "C", NA),
    loc5.2 = c("G", "G", NA),
    stringsAsFactors = FALSE
  )
  geno <- as_geno_data(df, id_col = "indiv.ID", group_col = NULL)
  res <- compare_replicates(geno, min_loci_compared = 1)

  # Check discrepancies table
  disc <- res$discrepancies
  expect_equal(nrow(disc), 3)
  expect_false("loc1" %in% disc$locus)
  expect_false("loc5" %in% disc$locus)

  # Check category for locus2 (allelic dropout)
  loc2_row <- disc[disc$locus == "loc2", ]
  expect_equal(loc2_row$discrepancy_type, "allelic_dropout")
  expect_equal(loc2_row$n_distinct_genotypes, 2)
  expect_equal(loc2_row$call_rep1, "A/T")
  expect_equal(loc2_row$call_rep2, "A/A")

  # Check category for locus3 (discordant)
  loc3_row <- disc[disc$locus == "loc3", ]
  expect_equal(loc3_row$discrepancy_type, "discordant")

  # Check category for locus4 (multiallelic_discordant)
  loc4_row <- disc[disc$locus == "loc4", ]
  expect_equal(loc4_row$discrepancy_type, "multiallelic_discordant")
})

test_that("compare_replicates computes ordination for n = 2 and n >= 3", {
  df <- data.frame(
    indiv.ID = c("s1", "s2", "s3", "s4"),
    l1.1 = c("A", "A", "G", "G"),
    l1.2 = c("A", "A", "G", "G"),
    l2.1 = c("C", "C", "C", "T"),
    l2.2 = c("C", "C", "C", "T"),
    stringsAsFactors = FALSE
  )
  geno <- as_geno_data(df, id_col = "indiv.ID", group_col = NULL)

  # n = 2
  res2 <- compare_replicates(geno, samples = c("s1", "s3"), min_loci_compared = 1)
  expect_equal(dim(res2$ordination$points), c(2, 2))
  expect_equal(res2$ordination$var_explained, c(100, 0))

  # n = 4
  res4 <- compare_replicates(geno, min_loci_compared = 1)
  expect_equal(dim(res4$ordination$points), c(4, 2))
  expect_true(all(!is.na(res4$ordination$points)))
  expect_true(sum(res4$ordination$var_explained) <= 100.01)
})

test_that("compare_replicates handles NA distances and warns when ordination is impossible", {
  geno <- synthetic_geno()
  # Request high min_loci_compared so some pairs become NA
  expect_warning(
    res <- compare_replicates(geno, samples = c("mom", "kid", "unrelated1"), min_loci_compared = 5),
    "Pairwise distance matrix contains NA values"
  )
  expect_null(res$ordination$points)
})

test_that("print.geno_replicates formats expected output", {
  geno <- synthetic_geno()
  res <- compare_replicates(geno, samples = c("kid", "dup_a", "dup_b"), min_loci_compared = 1)

  out <- capture.output(print(res))
  expect_true(any(grepl("<geno_replicates>", out)))
  expect_true(any(grepl("Pairwise concordance", out)))
  expect_true(any(grepl("Discrepant loci:", out)))
  expect_true(any(grepl("culprit loc", out)))
})

test_that("plot.geno_replicates works with base graphics and errors on missing ordination", {
  geno <- synthetic_geno()
  res <- compare_replicates(geno, samples = c("kid", "dup_a", "dup_b"), min_loci_compared = 1)

  # Base engine
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  pts <- plot(res, engine = "base")
  grDevices::dev.off()
  unlink(tf)

  expect_equal(dim(pts), c(3, 2))

  # plot_replicates alias
  tf2 <- tempfile(fileext = ".png")
  grDevices::png(tf2)
  pts2 <- plot_replicates(res, engine = "base")
  grDevices::dev.off()
  unlink(tf2)

  expect_equal(pts, pts2)

  # Error if ordination is NULL
  res_no_ord <- res
  res_no_ord$ordination$points <- NULL
  expect_error(plot(res_no_ord), "No ordination coordinates available")
})

test_that("plot.geno_replicates works with ggplot2 when installed", {
  skip_if_not_installed("ggplot2")
  geno <- synthetic_geno()
  res <- compare_replicates(geno, samples = c("kid", "dup_a", "dup_b"), min_loci_compared = 1)

  p <- plot(res, engine = "ggplot2")
  expect_s3_class(p, "ggplot")
})

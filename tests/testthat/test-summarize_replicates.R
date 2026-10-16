test_that("summarize_replicates input validation works", {
  geno <- synthetic_geno()
  set1 <- compare_replicates(geno, samples = c("kid", "dup_a"), min_loci_compared = 1)

  # Error if empty
  expect_error(summarize_replicates(), "Need at least one")

  # Error if non-geno_replicates passed
  expect_error(summarize_replicates(list(set1, "not_reps")), "must be `geno_replicates`")

  # Error if invalid thresholds
  expect_error(summarize_replicates(list(set1), min_concordance = 150), "between 0 and 100")
  expect_error(summarize_replicates(list(set1), locus_error_threshold = 2), "between 0 and 1")
})

test_that("summarize_replicates accepts both list and dot arguments", {
  geno <- synthetic_geno()
  set1 <- compare_replicates(geno, samples = c("kid", "dup_a"), min_loci_compared = 1)
  set2 <- compare_replicates(geno, samples = c("kid", "dup_b"), min_loci_compared = 1)

  # List input with names
  res1 <- summarize_replicates(list(FishA = set1, FishB = set2))
  expect_equal(res1$set_summary$set_name, c("FishA", "FishB"))

  # Dot arguments without names
  res2 <- summarize_replicates(set1, set2)
  expect_equal(res2$set_summary$set_name, c("set_1", "set_2"))
  expect_equal(res1$global$n_sets, 2)
})

test_that("summarize_replicates flags outlier sets and chronic culprit loci", {
  # Build 3 sets:
  # Set 1: kid vs dup_a (100% match)
  # Set 2: kid vs dup_b (80% match, mismatch at locus5)
  # Set 3: kid vs dup_b (80% match, mismatch at locus5)
  geno <- synthetic_geno()
  set1 <- compare_replicates(geno, samples = c("kid", "dup_a"), min_loci_compared = 1)
  set2 <- compare_replicates(geno, samples = c("kid", "dup_b"), min_loci_compared = 1)
  set3 <- compare_replicates(geno, samples = c("kid", "dup_b"), min_loci_compared = 1)

  trends <- summarize_replicates(
    list(Good = set1, Low1 = set2, Low2 = set3),
    min_concordance = 95,
    locus_error_threshold = 0.50
  )

  # Set summaries & outliers
  ss <- trends$set_summary
  expect_equal(nrow(ss), 3)
  expect_false(ss$flagged_outlier[ss$set_name == "Good"])
  expect_true(ss$flagged_outlier[ss$set_name == "Low1"])
  expect_true(ss$flagged_outlier[ss$set_name == "Low2"])

  # Locus trends
  lt <- trends$locus_trends
  # locus5 is discrepant in 2 of 3 sets (66.7%) -> flagged >= 50%
  loc5 <- lt[lt$locus == "locus5", ]
  expect_equal(loc5$n_sets_discrepant, 2)
  expect_equal(loc5$n_sets_evaluated, 3)
  expect_true(loc5$flagged)

  # Other loci clean
  clean_locs <- lt[lt$locus != "locus5", ]
  expect_true(all(!clean_locs$flagged))
  expect_true(all(clean_locs$discrepancy_rate == 0))

  # Global metrics
  expect_equal(trends$global$n_sets, 3)
  expect_equal(trends$global$n_sets_flagged, 2)
  expect_equal(trends$global$n_loci_flagged, 1)
  expect_equal(trends$global$n_loci_clean, 4)
  expect_equal(trends$global$pct_loci_clean, 80)
})

test_that("all_discrepancies combines records across sets", {
  geno <- synthetic_geno()
  set1 <- compare_replicates(geno, samples = c("kid", "dup_a"), min_loci_compared = 1)
  set2 <- compare_replicates(geno, samples = c("kid", "dup_b"), min_loci_compared = 1)

  trends <- summarize_replicates(list(SetA = set1, SetB = set2))
  ad <- trends$all_discrepancies

  expect_equal(nrow(ad), 1)
  expect_equal(ad$replicate_set, "SetB")
  expect_equal(ad$locus, "locus5")
})

test_that("print.replicate_trends outputs structured summary", {
  geno <- synthetic_geno()
  set1 <- compare_replicates(geno, samples = c("kid", "dup_a"), min_loci_compared = 1)
  set2 <- compare_replicates(geno, samples = c("kid", "dup_b"), min_loci_compared = 1)

  trends <- summarize_replicates(list(SetA = set1, SetB = set2), min_concordance = 95)
  out <- capture.output(print(trends))

  expect_true(any(grepl("<replicate_trends>", out)))
  expect_true(any(grepl("Global concordance", out)))
  expect_true(any(grepl("Panel robustness", out)))
  expect_true(any(grepl("flagged below 95% concordance", out)))
  expect_true(any(grepl("Chronic culprit loci", out)))
})

test_that("plot.replicate_trends works for both culprits and concordance", {
  geno <- synthetic_geno()
  set1 <- compare_replicates(geno, samples = c("kid", "dup_a"), min_loci_compared = 1)
  set2 <- compare_replicates(geno, samples = c("kid", "dup_b"), min_loci_compared = 1)
  trends <- summarize_replicates(list(SetA = set1, SetB = set2))

  # Culprits plot (base)
  tf1 <- tempfile(fileext = ".png")
  grDevices::png(tf1)
  p1 <- plot(trends, type = "culprits", engine = "base")
  grDevices::dev.off()
  unlink(tf1)
  expect_true(is.data.frame(p1))

  # Concordance plot (base)
  tf2 <- tempfile(fileext = ".png")
  grDevices::png(tf2)
  p2 <- plot(trends, type = "concordance", engine = "base")
  grDevices::dev.off()
  unlink(tf2)
  expect_true(is.data.frame(p2))

  # Alias plot_replicate_trends
  tf3 <- tempfile(fileext = ".png")
  grDevices::png(tf3)
  p3 <- plot_replicate_trends(trends, type = "culprits", engine = "base")
  grDevices::dev.off()
  unlink(tf3)
  expect_equal(p1, p3)
})

test_that("plot.replicate_trends works with ggplot2 when installed", {
  skip_if_not_installed("ggplot2")
  geno <- synthetic_geno()
  set1 <- compare_replicates(geno, samples = c("kid", "dup_a"), min_loci_compared = 1)
  set2 <- compare_replicates(geno, samples = c("kid", "dup_b"), min_loci_compared = 1)
  trends <- summarize_replicates(list(SetA = set1, SetB = set2))

  p1 <- plot(trends, type = "culprits", engine = "ggplot2")
  expect_s3_class(p1, "ggplot")

  p2 <- plot(trends, type = "concordance", engine = "ggplot2")
  expect_s3_class(p2, "ggplot")
})

test_that(".pairwise_counts computes hand-worked match/compared counts", {
  # locusA: always called; locusB: s3 missing (asymmetric missingness);
  # locusC: missing for every sample (zero-level edge case -- contributes
  # no columns to the sparse one-hot matrix at all).
  df <- data.frame(
    locusA.1 = c("A", "A", "G", "A"),
    locusA.2 = c("A", "A", "G", "G"),
    locusB.1 = c("C", "C", NA,  "C"),
    locusB.2 = c("C", "T", NA,  "C"),
    locusC.1 = c(NA, NA, NA, NA),
    locusC.2 = c(NA, NA, NA, NA),
    stringsAsFactors = FALSE
  )

  res <- .pairwise_counts(df, c("locusA", "locusB", "locusC"))

  expected_match <- matrix(c(
    2, 1, 0, 1,
    1, 2, 0, 0,
    0, 0, 1, 0,
    1, 0, 0, 2
  ), nrow = 4, byrow = TRUE)

  expected_compared <- matrix(c(
    2, 2, 1, 2,
    2, 2, 1, 2,
    1, 1, 1, 1,
    2, 2, 1, 2
  ), nrow = 4, byrow = TRUE)

  expect_equal(unname(res$match_count), expected_match)
  expect_equal(unname(res$compared_count), expected_compared)
})

test_that(".pairwise_counts handles zero loci and zero-level loci", {
  df <- data.frame(
    locusA.1 = c(NA_character_, NA_character_, NA_character_),
    locusA.2 = c(NA_character_, NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )

  res_zero_level <- .pairwise_counts(df, "locusA")
  expect_equal(unname(res_zero_level$match_count), matrix(0L, 3, 3))
  expect_equal(unname(res_zero_level$compared_count), matrix(0L, 3, 3))

  res_no_loci <- .pairwise_counts(df, character(0))
  expect_equal(unname(res_no_loci$match_count), matrix(0L, 3, 3))
  expect_equal(unname(res_no_loci$compared_count), matrix(0L, 3, 3))
})

# Reference copy of the pre-rewrite algorithm (O(n^2) per-locus outer()/
# matrix() comparisons), kept test-only as a fuzz-equivalence oracle for the
# sparse-crossprod .pairwise_counts() in R/utils.R -- guards against a future
# edit to this shared, non-obvious core helper silently changing behavior.
.pairwise_counts_reference_old <- function(x, loci) {
  n <- nrow(x)
  match_count <- matrix(0L, n, n)
  compared_count <- matrix(0L, n, n)

  for (locus in loci) {
    cols <- paste0(locus, c(".1", ".2"))
    key <- .canon_genotype_vec(x[[cols[1]]], x[[cols[2]]])
    present <- !is.na(key)

    both_present <- outer(present, present, "&")
    compared_count <- compared_count + both_present

    key_i <- matrix(key, n, n)
    key_j <- matrix(key, n, n, byrow = TRUE)
    match_count <- match_count + (both_present & (key_i == key_j))
  }

  list(match_count = match_count, compared_count = compared_count)
}

.random_geno_df <- function(n, L, n_alleles, missing_rate) {
  alleles <- sample(c("A", "C", "G", "T", "AC", "GT", "AAC", "GGT"), n_alleles)
  cols <- vector("list", L * 2)
  nm <- character(L * 2)
  for (li in seq_len(L)) {
    a1 <- sample(alleles, n, replace = TRUE)
    a2 <- sample(alleles, n, replace = TRUE)
    miss <- runif(n) < missing_rate
    a1[miss] <- NA
    a2[miss] <- NA
    cols[[2 * li - 1]] <- a1
    cols[[2 * li]] <- a2
    nm[2 * li - 1] <- paste0("locus", li, ".1")
    nm[2 * li] <- paste0("locus", li, ".2")
  }
  df <- as.data.frame(cols, stringsAsFactors = FALSE)
  names(df) <- nm
  df
}

test_that(".pairwise_counts matches the pre-rewrite reference algorithm across random shapes", {
  set.seed(20260811)
  configs <- list(
    list(n = 2,  L = 1,  n_alleles = 2, missing_rate = 0),
    list(n = 5,  L = 3,  n_alleles = 2, missing_rate = 0.3),
    list(n = 12, L = 8,  n_alleles = 4, missing_rate = 0.1),
    list(n = 25, L = 15, n_alleles = 6, missing_rate = 0.05),
    list(n = 40, L = 10, n_alleles = 3, missing_rate = 0.5),
    list(n = 15, L = 6,  n_alleles = 8, missing_rate = 0)
  )

  for (cfg in configs) {
    loci <- paste0("locus", seq_len(cfg$L))
    df <- .random_geno_df(cfg$n, cfg$L, cfg$n_alleles, cfg$missing_rate)

    old <- .pairwise_counts_reference_old(df, loci)
    new <- .pairwise_counts(df, loci)

    expect_identical(unname(new$match_count), unname(old$match_count))
    expect_identical(unname(new$compared_count), unname(old$compared_count))
  }
})

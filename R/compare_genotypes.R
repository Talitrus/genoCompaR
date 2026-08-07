.compare_loci <- function(x, row1, row2, loci, mendelian) {
  n <- length(loci)

  # Pull both samples' alleles for every locus in four vectorized reads
  # rather than one tibble subset per locus, then fill preallocated result
  # vectors and assemble a single tibble at the end. Building (and rbind-ing)
  # one small tibble per locus dominated runtime otherwise.
  a1 <- .allele_row(x, row1, loci, 1)
  a2 <- .allele_row(x, row1, loci, 2)
  b1 <- .allele_row(x, row2, loci, 1)
  b2 <- .allele_row(x, row2, loci, 2)

  status <- character(n)
  incompatible <- rep(NA, n)
  a_not_in_b <- rep(NA_character_, n)
  b_not_in_a <- rep(NA_character_, n)

  missing <- is.na(a1) | is.na(a2) | is.na(b1) | is.na(b2)
  # unordered genotype equality, vectorized
  equal <- (a1 == b1 & a2 == b2) | (a1 == b2 & a2 == b1)
  equal[missing] <- FALSE

  status[missing] <- "missing"
  status[!missing & equal] <- "match"
  status[!missing & !equal] <- "mismatch"

  if (mendelian) {
    incompatible[!missing] <- FALSE
    # only genuinely mismatching loci can be incompatible or have
    # unaccounted-for alleles, so the per-locus work is limited to those
    for (i in which(!missing & !equal)) {
      a <- c(a1[i], a2[i])
      b <- c(b1[i], b2[i])
      incompatible[i] <- length(intersect(a, b)) == 0
      d1 <- .multiset_setdiff(a, b)
      d2 <- .multiset_setdiff(b, a)
      if (length(d1) > 0) a_not_in_b[i] <- paste(d1, collapse = "/")
      if (length(d2) > 0) b_not_in_a[i] <- paste(d2, collapse = "/")
    }
  }

  tibble::tibble(
    locus = loci,
    sample1_a1 = a1, sample1_a2 = a2,
    sample2_a1 = b1, sample2_a2 = b2,
    status = status,
    mendelian_incompatible = incompatible,
    sample1_alleles_not_in_sample2 = a_not_in_b,
    sample2_alleles_not_in_sample1 = b_not_in_a
  )
}

#' Compare genotypes between two samples
#'
#' Performs a locus-by-locus comparison between two samples in a `geno_data`
#' object: how many loci both samples were successfully genotyped at ("loci
#' compared"), how many were only genotyped in one of the two ("unique
#' loci"), the number and percentage of loci with matching genotype calls,
#' and, optionally, a bi-directional Mendelian incompatibility screen useful
#' for checking a putative parent-offspring relationship.
#'
#' @section Mendelian incompatibility:
#' At each shared, non-missing locus, a genotype pair is *Mendelian
#' incompatible* when the two genotypes share no allele in common — i.e.
#' neither of sample 1's alleles appears in sample 2's genotype, and vice
#' versa. Allele-sharing is symmetric, so the incompatibility test itself
#' does not depend on which sample is treated as the putative parent.
#' "Bi-directional" refers to the *reporting*: for every incompatible locus,
#' [locus_detail()] records both which of sample 1's allele(s) are absent
#' from sample 2, and which of sample 2's allele(s) are absent from sample
#' 1, so the mismatch is visible from both samples' perspective. A single
#' incompatible locus is common due to genotyping error or mutation;
#' multiple incompatible loci across the panel is strong evidence against a
#' parent-offspring relationship.
#'
#' @param x A `geno_data` object, e.g. from [read_genotypes()].
#' @param sample1,sample2 Sample IDs (as found in `x`'s ID column) to
#'   compare.
#' @param mendelian Logical; if `TRUE` (the default) also compute
#'   bi-directional Mendelian incompatibility counts, appropriate for
#'   parentage screening. Set to `FALSE` for a plain identity/match-rate
#'   comparison (e.g. duplicate/twin checking) without the parentage
#'   framing.
#' @param output One of `"pretty"` (default) or `"tibble"`. `"pretty"`
#'   prints a human-readable report to the console and returns the summary
#'   tibble invisibly. `"tibble"` returns the one-row summary tibble
#'   directly with no console output — convenient when looping over many
#'   pairs.
#'
#' @return A `geno_comparison` object (subclass of tibble) with one summary
#'   row and columns `sample1`, `sample2`, `n_loci_total`, `n_loci_compared`,
#'   `n_loci_unique` (genotyped in only one sample), `n_match`, `n_mismatch`,
#'   `pct_match`, and, when `mendelian = TRUE`, `n_mendelian_incompatible`
#'   and `pct_mendelian_incompatible`. The full per-locus comparison table is
#'   attached and retrievable with [locus_detail()].
#'
#' @seealso [find_duplicate_genotypes()] to scan an entire dataset for
#'   matching pairs; [locus_detail()] to pull the full per-locus table.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' compare_genotypes(geno, "CH-00001_S1", "CH-00002_S2")
#'
#' result <- compare_genotypes(
#'   geno, "CH-00001_S1", "CH-00002_S2",
#'   output = "tibble"
#' )
#' result
#' locus_detail(result)
#'
#' @export
compare_genotypes <- function(x,
                               sample1,
                               sample2,
                               mendelian = TRUE,
                               output = c("pretty", "tibble")) {
  output <- match.arg(output)
  if (!inherits(x, "geno_data")) {
    stop("`x` must be a geno_data object (see read_genotypes()).", call. = FALSE)
  }
  if (identical(sample1, sample2)) {
    stop("`sample1` and `sample2` must be different samples.", call. = FALSE)
  }

  row1 <- .row_index(x, sample1)
  row2 <- .row_index(x, sample2)
  loci <- get_loci(x)

  detail <- .compare_loci(x, row1, row2, loci, mendelian = mendelian)

  n_total    <- length(loci)
  n_compared <- sum(detail$status != "missing")
  n_unique   <- sum(detail$status == "missing")
  n_match    <- sum(detail$status == "match")
  n_mismatch <- sum(detail$status == "mismatch")
  pct_match  <- if (n_compared > 0) 100 * n_match / n_compared else NA_real_

  summary <- tibble::tibble(
    sample1 = sample1,
    sample2 = sample2,
    n_loci_total = n_total,
    n_loci_compared = n_compared,
    n_loci_unique = n_unique,
    n_match = n_match,
    n_mismatch = n_mismatch,
    pct_match = pct_match
  )

  if (mendelian) {
    n_incompat <- sum(detail$mendelian_incompatible %in% TRUE)
    summary$n_mendelian_incompatible <- n_incompat
    summary$pct_mendelian_incompatible <-
      if (n_compared > 0) 100 * n_incompat / n_compared else NA_real_
  }

  result <- structure(
    summary,
    class = c("geno_comparison", class(summary)),
    locus_detail = detail,
    mendelian = mendelian
  )

  if (output == "pretty") {
    print(result)
    return(invisible(result))
  }
  .as_plain_summary(result)
}

#' Retrieve the per-locus detail table from a comparison result
#'
#' @param x A `geno_comparison` object returned by [compare_genotypes()].
#' @return A tibble with one row per locus, giving both samples' alleles,
#'   match/mismatch/missing `status`, and (when the comparison was run with
#'   `mendelian = TRUE`) Mendelian incompatibility flags and the alleles
#'   unaccounted for in each direction.
#' @export
locus_detail <- function(x) {
  d <- attr(x, "locus_detail")
  if (is.null(d)) {
    stop(
      "`x` has no locus_detail attribute; is it a geno_comparison object ",
      "from compare_genotypes()?",
      call. = FALSE
    )
  }
  d
}

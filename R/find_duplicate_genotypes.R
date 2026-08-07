#' Scan a dataset for matching / duplicate genotypes
#'
#' Compares every pair of samples in a `geno_data` object and flags pairs
#' whose genotypes match at or above `threshold` of their shared
#' (non-missing) loci — i.e. potential duplicate samples, technical
#' replicates, or clones ("matchy pairs"). Useful for spotting sample
#' mix-ups, re-sampled individuals, or accidental re-runs before analysis.
#'
#' @param x A `geno_data` object, e.g. from [read_genotypes()].
#' @param threshold Minimum proportion (0-1) of compared loci that must
#'   match for a pair to be reported. Defaults to `0.9`.
#' @param min_loci_compared Minimum number of loci both samples must have
#'   been genotyped at for a pair to be considered; guards against small
#'   numbers of compared loci producing spuriously high match percentages.
#'   Defaults to `30`.
#' @param loci Character vector of loci to use. Defaults to all loci in `x`.
#'
#' @return A `geno_duplicates` object (subclass of tibble), one row per
#'   flagged pair, sorted by descending `pct_match`, with columns
#'   `sample1`, `sample2`, `n_loci_compared`, `n_match`, `n_mismatch`,
#'   `pct_match` (and `group1`/`group2` if `x` has a group column).
#'
#' @details
#' This runs in O(n^2 x L) time/memory in the number of samples (n) and
#' loci (L), using vectorized matrix comparisons rather than looping over
#' individual pairs. For very large datasets (many thousands of samples),
#' consider filtering `x` down (e.g. to one group/run at a time, or a subset
#' of high-confidence loci via `loci`) before calling this function.
#'
#' @seealso [compare_genotypes()] for a detailed look at one pair.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' dupes <- find_duplicate_genotypes(geno, threshold = 0.95)
#' dupes
#'
#' @export
find_duplicate_genotypes <- function(x,
                                      threshold = 0.9,
                                      min_loci_compared = 30,
                                      loci = get_loci(x)) {
  if (!inherits(x, "geno_data")) {
    stop("`x` must be a geno_data object (see read_genotypes()).", call. = FALSE)
  }
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("`threshold` must be a number between 0 and 1.", call. = FALSE)
  }
  n <- nrow(x)
  if (n < 2) {
    stop("Need at least 2 samples to compare.", call. = FALSE)
  }

  counts <- .pairwise_counts(x, loci)
  match_count <- counts$match_count
  compared_count <- counts$compared_count

  pct_match <- ifelse(compared_count > 0, match_count / compared_count, NA_real_)

  ids <- x[[.id_col(x)]]
  pair_idx <- which(upper.tri(pct_match), arr.ind = TRUE)

  keep <- !is.na(pct_match[pair_idx]) &
    pct_match[pair_idx] >= threshold &
    compared_count[pair_idx] >= min_loci_compared
  pair_idx <- pair_idx[keep, , drop = FALSE]

  result <- tibble::tibble(
    sample1 = ids[pair_idx[, 1]],
    sample2 = ids[pair_idx[, 2]],
    n_loci_compared = compared_count[pair_idx],
    n_match = match_count[pair_idx],
    n_mismatch = compared_count[pair_idx] - match_count[pair_idx],
    pct_match = 100 * pct_match[pair_idx]
  )

  gc <- .group_col(x)
  if (!is.null(gc)) {
    result$group1 <- x[[gc]][pair_idx[, 1]]
    result$group2 <- x[[gc]][pair_idx[, 2]]
  }

  result <- result[order(-result$pct_match), ]

  structure(
    result,
    class = c("geno_duplicates", class(result)),
    n_samples = n,
    n_loci = length(loci),
    threshold = threshold
  )
}

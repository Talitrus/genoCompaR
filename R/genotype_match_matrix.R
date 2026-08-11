#' Pairwise genotype match-rate matrix
#'
#' Computes the percent-matching genotype rate between every pair of
#' samples in a `geno_data` object, returning a symmetric sample x sample
#' matrix. Useful as input to clustering (e.g. `stats::hclust()` on a
#' distance derived from `1 - matrix / 100`) or a heatmap, or for eyeballing
#' overall relatedness/duplication structure rather than only the pairs
#' that clear a fixed threshold (see [find_duplicate_genotypes()] for that).
#'
#' @param x A `geno_data` object, e.g. from [read_genotypes()].
#' @param min_loci_compared Minimum number of loci both samples in a pair
#'   must have been genotyped at for that cell to be computed; pairs below
#'   this are set to `NA` rather than reporting a match rate based on very
#'   little data. Defaults to `30`, matching
#'   [find_duplicate_genotypes()].
#' @param loci Character vector of loci to use. Defaults to all loci in `x`.
#'
#' @return A numeric matrix of class `geno_match_matrix`, `n x n` (n =
#'   `nrow(x)`), row/column-named by the sample ID column, with entries
#'   equal to the percentage of compared loci that match between each pair
#'   of samples. The diagonal is `NA` (self-comparison is not meaningful).
#'   The underlying compared-loci counts are attached as attribute
#'   `"n_compared"`.
#'
#' @details
#' Like [find_duplicate_genotypes()], this runs in O(n^2 x L) time/memory,
#' via a sparse matrix multiply rather than per-locus O(n^2) comparisons;
#' this keeps it fast in practice (thousands of samples in seconds, not
#' minutes). Filtering `x` first (e.g. to one group/run, or a subset of
#' loci) remains an option but isn't required for performance at the sizes
#' this package targets.
#'
#' Real panels routinely contain sample pairs that share too few typed loci
#' to score, so **the matrix normally contains `NA` cells** (in addition to
#' the always-`NA` diagonal). [stats::hclust()] rejects `NA` distances, so
#' clustering requires deciding what an unscorable pair means first - see
#' the example below.
#'
#' @seealso [find_duplicate_genotypes()] to get only pairs above a match
#'   threshold as a tidy table; [compare_genotypes()] for a detailed look
#'   at one pair.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' mm <- genotype_match_matrix(geno)
#' mm[1:5, 1:5]
#'
#' # Hierarchical clustering on genotype similarity. Pairs sharing too few
#' # typed loci come back NA, which hclust() cannot accept; treating them as
#' # maximally distant is the conservative reading (no evidence they match).
#' d <- stats::as.dist(1 - mm / 100)
#' d[is.na(d)] <- 1
#' plot(stats::hclust(d), labels = FALSE)
#'
#' @export
genotype_match_matrix <- function(x,
                                   min_loci_compared = 30,
                                   loci = get_loci(x)) {
  if (!inherits(x, "geno_data")) {
    stop("`x` must be a geno_data object (see read_genotypes()).", call. = FALSE)
  }
  n <- nrow(x)
  if (n < 2) {
    stop("Need at least 2 samples to compare.", call. = FALSE)
  }

  counts <- .pairwise_counts(x, loci)
  ok <- counts$compared_count > 0 & counts$compared_count >= min_loci_compared
  pct_match <- ifelse(ok, 100 * counts$match_count / counts$compared_count, NA_real_)
  diag(pct_match) <- NA_real_

  ids <- x[[.id_col(x)]]
  dimnames(pct_match) <- list(ids, ids)

  structure(
    pct_match,
    class = c("geno_match_matrix", class(pct_match)),
    n_compared = counts$compared_count
  )
}

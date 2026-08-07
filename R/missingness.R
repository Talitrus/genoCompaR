#' Per-locus missing-data summary
#'
#' Counts, for each locus, how many samples have a missing call (either
#' allele column `NA`).
#'
#' @param x A `geno_data` object, e.g. from [read_genotypes()].
#' @return A tibble with one row per locus: `locus`, `n_missing`, `n_typed`,
#'   `pct_missing`.
#' @seealso [sample_missingness()]
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' locus_missingness(geno)
#'
#' @export
locus_missingness <- function(x) {
  if (!inherits(x, "geno_data")) {
    stop("`x` must be a geno_data object (see read_genotypes()).", call. = FALSE)
  }
  loci <- get_loci(x)
  n <- nrow(x)

  n_missing <- vapply(loci, function(l) {
    cols <- paste0(l, c(".1", ".2"))
    sum(is.na(x[[cols[1]]]) | is.na(x[[cols[2]]]))
  }, integer(1), USE.NAMES = FALSE)

  tibble::tibble(
    locus = loci,
    n_missing = n_missing,
    n_typed = n - n_missing,
    pct_missing = 100 * n_missing / n
  )
}

#' Per-sample missing-data summary
#'
#' Counts, for each sample, how many loci have a missing call.
#'
#' @param x A `geno_data` object, e.g. from [read_genotypes()].
#' @return A tibble with one row per sample: the ID column, `n_loci_missing`,
#'   `n_loci_typed`, `pct_missing`.
#' @seealso [locus_missingness()]
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' sample_missingness(geno)
#'
#' @export
sample_missingness <- function(x) {
  if (!inherits(x, "geno_data")) {
    stop("`x` must be a geno_data object (see read_genotypes()).", call. = FALSE)
  }
  loci <- get_loci(x)
  n_loci <- length(loci)

  a1 <- as.matrix(x[paste0(loci, ".1")])
  a2 <- as.matrix(x[paste0(loci, ".2")])
  # as.integer to match locus_missingness()'s integer counts
  n_missing <- as.integer(rowSums(is.na(a1) | is.na(a2)))

  id_name <- .id_col(x)
  out <- data.frame(
    id = x[[id_name]],
    n_loci_missing = n_missing,
    n_loci_typed = n_loci - n_missing,
    pct_missing = 100 * n_missing / n_loci,
    stringsAsFactors = FALSE
  )
  names(out)[1] <- id_name
  tibble::as_tibble(out)
}

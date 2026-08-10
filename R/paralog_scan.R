#' Scan loci for heterozygote excess (collapsed paralogs)
#'
#' Compares each locus's observed heterozygosity against the Hardy-Weinberg
#' expectation computed from that locus's *own* allele frequencies, and flags
#' loci that are called heterozygous far more often than any frequency
#' distribution permits. That is the signature of two paralogous loci being
#' genotyped as one: every individual carries one copy of each paralog, so
#' nearly every individual is called heterozygous.
#'
#' @param x A `geno_data` object (e.g. from [read_genotypes()]), a
#'   `geno_conformed` table from [conform_genotypes()], or any wide data frame
#'   with two allele columns per locus (`<locus>.1`, `<locus>.2`).
#' @param min_n Loci with fewer than this many complete genotypes are skipped;
#'   the statistic is meaningless on small samples. Defaults to `50`.
#' @param threshold Excess (`obs_het - exp_het`) at or above which a locus is
#'   flagged. Defaults to `0.25`.
#' @param loci Character vector of loci to scan. Defaults to every locus in
#'   `x`.
#'
#' @return A `geno_paralog_scan` object (subclass of tibble), one row per
#'   locus with enough data, sorted by descending `excess`, with columns
#'   `locus`, `n`, `alleles`, `obs_het`, `exp_het`, `excess`, `top_alleles`
#'   (the three commonest alleles and their counts), and `flagged`.
#'
#' @details
#' Testing against the locus's own expectation rather than a flat 0.5 is what
#' keeps legitimately diverse multi-allelic microhaplotypes from being
#' condemned: a biallelic locus caps at 0.5 expected heterozygosity, but a
#' microhap with six alleles can sit well above it honestly. Only the gap
#' between observed and expected is evidence of a problem.
#'
#' Nothing else in the package's QC sees this fault. The individual calls are
#' well-formed haplotypes of the right width, so both the width check and the
#' VCF allele check in [conform_genotypes()] pass them.
#'
#' @seealso [conform_genotypes()], which runs this scan as part of its QC
#'   report; [locus_missingness()] for the other per-locus summary.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' scan <- paralog_scan(geno, min_n = 30)
#' scan
#'
#' # the flagged loci alone
#' scan[scan$flagged, c("locus", "obs_het", "exp_het", "excess")]
#'
#' @export
paralog_scan <- function(x, min_n = 50L, threshold = 0.25, loci = NULL) {
  if (!is.data.frame(x)) {
    stop(
      "`x` must be a data frame, a geno_data object (see read_genotypes()), ",
      "or a geno_conformed table (see conform_genotypes()).",
      call. = FALSE
    )
  }
  if (is.null(loci)) {
    loci <- if (inherits(x, "geno_data")) {
      get_loci(x)
    } else {
      .split_genotype_columns(names(x))$loci
    }
  }
  if (length(loci) == 0) {
    stop("No locus allele-column pairs ('<locus>.1' / '<locus>.2') found.", call. = FALSE)
  }

  rows <- lapply(loci, function(locus) {
    cols <- paste0(locus, c(".1", ".2"))
    if (!all(cols %in% names(x))) return(NULL)
    a <- x[[cols[1]]]
    b <- x[[cols[2]]]
    ok <- .is_called(a) & .is_called(b)
    n <- sum(ok)
    if (n < min_n) return(NULL)

    freq <- table(c(a[ok], b[ok]))
    p <- freq / sum(freq)
    top <- utils::head(sort(freq, decreasing = TRUE), 3L)
    data.frame(
      locus = locus,
      n = n,
      alleles = length(freq),
      obs_het = sum(a[ok] != b[ok]) / n,
      exp_het = 1 - sum(p^2),
      top_alleles = paste(sprintf("%s=%d", names(top), as.integer(top)),
                          collapse = " "),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(
      locus = character(0), n = integer(0), alleles = integer(0),
      obs_het = numeric(0), exp_het = numeric(0), top_alleles = character(0),
      stringsAsFactors = FALSE
    )
  }
  out$excess <- out$obs_het - out$exp_het
  out$flagged <- out$excess > threshold
  out <- out[order(-out$excess),
             c("locus", "n", "alleles", "obs_het", "exp_het", "excess",
               "top_alleles", "flagged")]

  result <- tibble::as_tibble(out)
  structure(
    result,
    class = c("geno_paralog_scan", class(result)),
    n_loci = length(loci),
    min_n = min_n,
    threshold = threshold
  )
}

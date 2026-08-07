.trio_loci <- function(x, off_row, p1_row, p2_row, loci) {
  # Vectorized for the same reason as .compare_loci(): whole-column reads and
  # a single tibble at the end, instead of one tibble per locus.
  o1 <- .allele_row(x, off_row, loci, 1)
  o2 <- .allele_row(x, off_row, loci, 2)
  p11 <- .allele_row(x, p1_row, loci, 1)
  p12 <- .allele_row(x, p1_row, loci, 2)
  p21 <- .allele_row(x, p2_row, loci, 1)
  p22 <- .allele_row(x, p2_row, loci, 2)

  missing <- is.na(o1) | is.na(o2) | is.na(p11) | is.na(p12) |
    is.na(p21) | is.na(p22)

  # "allele o can come from parent p" reduces to membership in p's two alleles
  o1_from_p1 <- o1 == p11 | o1 == p12
  o1_from_p2 <- o1 == p21 | o1 == p22
  o2_from_p1 <- o2 == p11 | o2 == p12
  o2_from_p2 <- o2 == p21 | o2 == p22

  # homozygous offspring: the single allele must be available from both
  # parents; heterozygous: some assignment of one allele to each must work
  homozygous <- !missing & o1 == o2
  compatible <- ifelse(
    homozygous,
    o1_from_p1 & o1_from_p2,
    (o1_from_p1 & o2_from_p2) | (o2_from_p1 & o1_from_p2)
  )

  status <- ifelse(missing, "missing",
                   ifelse(compatible, "compatible", "incompatible"))

  tibble::tibble(
    locus = loci,
    offspring_a1 = o1, offspring_a2 = o2,
    parent1_a1 = p11, parent1_a2 = p12,
    parent2_a1 = p21, parent2_a2 = p22,
    status = status
  )
}

#' Trio-based Mendelian parentage test
#'
#' Tests, locus by locus, whether an offspring's genotype is fully
#' explainable by inheriting one allele from `parent1` and one allele from
#' `parent2`. This is a stronger test than running two independent pairwise
#' [compare_genotypes()] duo checks: a duo check only requires the offspring
#' to share *some* allele with each candidate parent, which can pass even
#' when no consistent one-allele-from-each assignment actually exists across
#' both parents at once. `check_trio()` requires a valid assignment.
#'
#' @section How compatibility is determined:
#' At each locus where all three samples are genotyped, let O be the
#' offspring's two alleles, P1 parent1's two alleles, and P2 parent2's two
#' alleles.
#' * If the offspring is homozygous (both alleles equal to `o`), the locus
#'   is compatible iff `o` appears in both P1 and P2.
#' * If the offspring is heterozygous (alleles `o1` and `o2`), the locus is
#'   compatible iff there is a way to assign one offspring allele to
#'   parent1 and the other to parent2: (`o1` in P1 and `o2` in P2), or
#'   (`o2` in P1 and `o1` in P2).
#'
#' `parent1` and `parent2` are interchangeable — the test does not assume
#' which one is the mother or the father. A single incompatible locus is
#' common due to genotyping error or mutation; multiple incompatible loci
#' across the panel is strong evidence against this trio being a true
#' family.
#'
#' @param x A `geno_data` object, e.g. from [read_genotypes()].
#' @param offspring,parent1,parent2 Sample IDs (as found in `x`'s ID
#'   column). `parent1`/`parent2` order doesn't matter.
#' @param output One of `"pretty"` (default) or `"tibble"`, as in
#'   [compare_genotypes()].
#'
#' @return A `geno_trio` object (subclass of tibble) with one summary row:
#'   `offspring`, `parent1`, `parent2`, `n_loci_total`, `n_loci_compared`
#'   (all three samples genotyped), `n_loci_unique` (at least one sample
#'   missing at that locus), `n_compatible`, `n_incompatible`,
#'   `pct_compatible`. The full per-locus table is attached and retrievable
#'   with [locus_detail()].
#'
#' @seealso [compare_genotypes()] for a pairwise duo check.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' check_trio(geno, "CH-00001_S1", "CH-00002_S2", "CH-00003_S3")
#'
#' @export
check_trio <- function(x,
                        offspring,
                        parent1,
                        parent2,
                        output = c("pretty", "tibble")) {
  output <- match.arg(output)
  if (!inherits(x, "geno_data")) {
    stop("`x` must be a geno_data object (see read_genotypes()).", call. = FALSE)
  }
  if (anyDuplicated(c(offspring, parent1, parent2))) {
    stop(
      "`offspring`, `parent1`, and `parent2` must all be different samples.",
      call. = FALSE
    )
  }

  off_row <- .row_index(x, offspring)
  p1_row  <- .row_index(x, parent1)
  p2_row  <- .row_index(x, parent2)
  loci <- get_loci(x)

  detail <- .trio_loci(x, off_row, p1_row, p2_row, loci)

  n_total    <- length(loci)
  n_compared <- sum(detail$status != "missing")
  n_unique   <- sum(detail$status == "missing")
  n_compat   <- sum(detail$status == "compatible")
  n_incompat <- sum(detail$status == "incompatible")
  pct_compat <- if (n_compared > 0) 100 * n_compat / n_compared else NA_real_

  summary <- tibble::tibble(
    offspring = offspring,
    parent1 = parent1,
    parent2 = parent2,
    n_loci_total = n_total,
    n_loci_compared = n_compared,
    n_loci_unique = n_unique,
    n_compatible = n_compat,
    n_incompatible = n_incompat,
    pct_compatible = pct_compat
  )

  result <- structure(
    summary,
    class = c("geno_trio", class(summary)),
    locus_detail = detail
  )

  if (output == "pretty") {
    print(result)
    return(invisible(result))
  }
  .as_plain_summary(result)
}

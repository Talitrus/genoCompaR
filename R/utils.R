# Internal helpers shared across genoCompaR functions.
# Not exported except where noted.

.validate_locus_columns <- function(cols) {
  if (length(cols) == 0 || length(cols) %% 2 != 0) {
    stop(
      "Expected an even number of allele columns (each locus needs a '.1' ",
      "and '.2' column); found ", length(cols), ".",
      call. = FALSE
    )
  }
  bad <- !grepl("\\.(1|2)$", cols)
  if (any(bad)) {
    stop(
      "Locus columns must end in '.1' or '.2'; offending columns: ",
      paste(cols[bad], collapse = ", "),
      call. = FALSE
    )
  }
  base <- sub("\\.(1|2)$", "", cols)
  loci <- unique(base)
  for (l in loci) {
    pair <- paste0(l, c(".1", ".2"))
    if (!all(pair %in% cols)) {
      stop(
        "Locus '", l, "' is missing one of its allele columns (expected ",
        paste(pair, collapse = " and "), ").",
        call. = FALSE
      )
    }
  }
  loci
}

# Drops the report-printing class (geno_comparison / geno_trio) so the object
# prints and behaves as an ordinary tibble, while keeping the locus_detail
# attribute so locus_detail() still works. Used for output = "tibble": those
# classes exist only to drive the pretty console report, and leaving them on
# would make a returned tibble render as a one-pair report - actively
# misleading once several results are row-bound together.
.as_plain_summary <- function(x) {
  class(x) <- setdiff(class(x), c("geno_comparison", "geno_trio"))
  x
}

.geno_attr <- function(x, which) {
  v <- attr(x, which)
  if (is.null(v)) {
    stop("`x` must be a geno_data object (see read_genotypes()).", call. = FALSE)
  }
  v
}

#' Locus names carried by a `geno_data` object
#'
#' @param x A `geno_data` object, e.g. from [read_genotypes()].
#' @return Character vector of locus names (without the `.1`/`.2` suffix).
#' @export
get_loci <- function(x) {
  .geno_attr(x, "loci")
}

.id_col <- function(x) .geno_attr(x, "id_col")

.group_col <- function(x) {
  gc <- attr(x, "group_col")
  if (is.null(gc) || is.na(gc)) NULL else gc
}

.row_index <- function(x, sample_id) {
  ids <- x[[.id_col(x)]]
  idx <- match(sample_id, ids)
  if (is.na(idx)) {
    stop(
      "Sample '", sample_id, "' not found in column '", .id_col(x), "'.",
      call. = FALSE
    )
  }
  idx
}

.allele_pair <- function(x, row, locus) {
  cols <- paste0(locus, c(".1", ".2"))
  unname(unlist(x[row, cols]))
}

# One sample's allele `which` (1 or 2) across many loci, as a plain character
# vector. Reads whole columns rather than subsetting per locus, which is what
# makes the per-locus comparison builders fast.
.allele_row <- function(x, row, loci, which) {
  cols <- paste0(loci, ".", which)
  vapply(cols, function(cn) x[[cn]][row], character(1), USE.NAMES = FALSE)
}

# Unordered (multiset) equality of two length-2 allele vectors.
.multiset_equal <- function(a, b) {
  identical(sort(a), sort(b))
}

# Elements of `a` (with multiplicity) that have no matching partner in `b`.
.multiset_setdiff <- function(a, b) {
  pool <- b
  out <- character(0)
  for (v in a) {
    hit <- match(v, pool)
    if (is.na(hit)) {
      out <- c(out, v)
    } else {
      pool <- pool[-hit]
    }
  }
  out
}

# Vectorized canonical genotype key for a locus across all samples: sorted
# "allele1/allele2" string, or NA if either allele is missing. Used for the
# fast all-pairs scan shared by find_duplicate_genotypes() and
# genotype_match_matrix().
.canon_genotype_vec <- function(a1, a2) {
  swap <- a1 > a2
  swap[is.na(swap)] <- FALSE
  lo <- ifelse(swap, a2, a1)
  hi <- ifelse(swap, a1, a2)
  key <- paste(lo, hi, sep = "/")
  key[is.na(a1) | is.na(a2)] <- NA_character_
  key
}

# All-pairs match_count / compared_count matrices (n x n, unnamed) for a set
# of loci. O(n^2 x L) via vectorized matrix comparisons, one locus at a time.
.pairwise_counts <- function(x, loci) {
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

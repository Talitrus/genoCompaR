# Internal helpers shared across genoCompaR functions.
# Not exported except where noted.

# TRUE where a genotype call is actually present. Call sheets have used all
# three spellings of "not called" - a real NA, an empty cell, and the literal
# string "NA" - so anything reading a file with colClasses = "character" (which
# conform_genotypes() must, to keep numeric-looking haplotypes from becoming
# integers) has to treat them alike.
.is_called <- function(x) !is.na(x) & x != "" & x != "NA"

# Read a call sheet with every column left as character and no NA conversion.
# conform_genotypes() needs both: colClasses keeps all-numeric marker columns
# from becoming integers (nchar() would then report digit counts instead of
# haplotype widths), and leaving "" / "NA" untouched keeps the raw file's own
# spelling of "uncalled" visible to .is_called().
.read_chr <- function(path, ...) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                  colClasses = "character", ...)
}

# Split an arbitrary wide genotype table's column names into locus base names
# and everything else. Unlike .validate_locus_columns() this is tolerant: it
# takes whatever pairs up and treats the rest as metadata, which is what
# paralog_scan() needs on a data frame that has no geno_data attributes and
# whose ID/group columns are not known in advance.
.split_genotype_columns <- function(cols) {
  firsts <- cols[grepl("\\.1$", cols)]
  seconds <- cols[grepl("\\.2$", cols)]
  loci <- intersect(sub("\\.1$", "", firsts), sub("\\.2$", "", seconds))
  genotype_cols <- c(paste0(loci, ".1"), paste0(loci, ".2"))
  list(loci = sort(loci), meta = setdiff(cols, genotype_cols))
}

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
# of loci. Recodes each locus's canonical genotype key into integer levels
# scoped into one global column space, builds a single sparse n x (total
# levels) one-hot indicator matrix H, and gets match_count from one
# Matrix::tcrossprod(H) sparse matrix multiply instead of an O(n^2) R-level
# comparison per locus. compared_count uses a plain dense tcrossprod() of
# the n x L presence matrix -- presence is usually ~95% dense, so routing
# it through a sparse conversion is slower than base R's BLAS-backed
# tcrossprod(), not faster.
.pairwise_counts <- function(x, loci) {
  n <- nrow(x)
  L <- length(loci)

  i_list <- vector("list", L)
  j_list <- vector("list", L)
  presence <- matrix(0, n, L)

  col_offset <- 0L
  for (li in seq_len(L)) {
    locus <- loci[li]
    cols <- paste0(locus, c(".1", ".2"))
    key <- .canon_genotype_vec(x[[cols[1]]], x[[cols[2]]])
    present <- !is.na(key)
    presence[, li] <- present

    f <- factor(key[present])
    n_levels <- nlevels(f)

    i_list[[li]] <- which(present)
    j_list[[li]] <- col_offset + as.integer(f)
    col_offset <- col_offset + n_levels
  }

  ii <- unlist(i_list, use.names = FALSE)
  jj <- unlist(j_list, use.names = FALSE)
  H <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, col_offset))
  match_count <- as.matrix(Matrix::tcrossprod(H))
  storage.mode(match_count) <- "integer"

  compared_count <- tcrossprod(presence)
  storage.mode(compared_count) <- "integer"

  list(match_count = match_count, compared_count = compared_count)
}

#' Build a `geno_data` object from a data frame
#'
#' Validates and classes an existing data frame / tibble that is already in
#' wide genotype layout (one row per sample, two allele columns per locus
#' named `<locus>.1` and `<locus>.2`). Use this when the genotypes did not
#' come straight from a single CSV read by [read_genotypes()] - for example
#' after subsetting rows, selecting a locus panel, or assembling a table
#' from another source.
#'
#' @param df A data frame or tibble in wide genotype layout.
#' @param id_col Name of the column holding unique sample identifiers.
#'   Defaults to `"indiv.ID"`.
#' @param group_col Name of an optional grouping column, or `NULL` for none.
#'   Defaults to `"group"` when such a column is present.
#' @param loci Character vector of locus names to keep. Defaults to every
#'   locus found in `df`.
#'
#' @return A `geno_data` object (see [read_genotypes()]).
#'
#' @seealso [read_genotypes()] to read from a file; [merge_genotypes()] to
#'   combine two genotype sets.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#'
#' # keep three samples and the first 20 loci
#' small <- as_geno_data(
#'   as.data.frame(geno)[1:3, ],
#'   loci = get_loci(geno)[1:20]
#' )
#' small
#'
#' @export
as_geno_data <- function(df,
                          id_col = "indiv.ID",
                          group_col = "group",
                          loci = NULL) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame or tibble.", call. = FALSE)
  }
  if (!id_col %in% names(df)) {
    stop("id_col '", id_col, "' not found in `df`.", call. = FALSE)
  }
  ids <- as.character(df[[id_col]])
  if (anyDuplicated(ids)) {
    dupes <- unique(ids[duplicated(ids)])
    stop(
      "Duplicate sample IDs in id_col '", id_col, "': ",
      paste(utils::head(dupes, 5), collapse = ", "),
      if (length(dupes) > 5) ", ..." else "",
      call. = FALSE
    )
  }

  has_group <- !is.null(group_col) && group_col %in% names(df)
  meta_cols <- c(id_col, if (has_group) group_col)

  if (is.null(loci)) {
    loci <- .validate_locus_columns(setdiff(names(df), meta_cols))
  } else {
    wanted <- as.vector(rbind(paste0(loci, ".1"), paste0(loci, ".2")))
    missing_cols <- setdiff(wanted, names(df))
    if (length(missing_cols) > 0) {
      stop(
        "`df` is missing allele column(s): ",
        paste(utils::head(missing_cols, 5), collapse = ", "),
        if (length(missing_cols) > 5) ", ..." else "",
        call. = FALSE
      )
    }
  }

  locus_cols <- as.vector(rbind(paste0(loci, ".1"), paste0(loci, ".2")))
  out <- tibble::as_tibble(df[, c(meta_cols, locus_cols), drop = FALSE])
  out[[id_col]] <- ids

  structure(
    out,
    class = c("geno_data", class(out)),
    id_col = id_col,
    group_col = if (has_group) group_col else NA_character_,
    loci = loci
  )
}

#' Combine two genotype sets on their shared loci
#'
#' Row-binds two `geno_data` objects, restricting to the loci they have in
#' common. This is the usual way to compare samples that arrived in separate
#' files - for example a broodstock reference panel and a set of unknown
#' individuals, which frequently use different (though overlapping) locus
#' panels and different sample-ID column names.
#'
#' @param x,y `geno_data` objects, e.g. from [read_genotypes()] or
#'   [as_geno_data()].
#' @param id_col Name for the sample ID column of the result. Defaults to
#'   `x`'s ID column; `y`'s ID column is renamed to match.
#' @param loci Character vector of loci to keep. Defaults to the shared loci,
#'   or every locus in either input when `keep_unmatched_loci = TRUE`.
#' @param keep_unmatched_loci Keep loci present in only one input as extra
#'   columns, NA-filled for the rows of the input that lacks them, instead of
#'   restricting to the shared panel. `FALSE` by default: for the duo/trio
#'   comparisons this function normally feeds, a locus one side never
#'   assayed can never be compared across the merge, so keeping it just
#'   inflates `n_loci_total` and dilutes per-sample call rates computed over
#'   the merged table.
#'
#' @return A `geno_data` object containing every row of `x` then every row
#'   of `y`, over the requested loci. The group column is retained only when
#'   both inputs have one.
#'
#' @details
#' Sample IDs must be unique across the combined result, since the rest of
#' the package assumes one row per sample; an ID appearing in both inputs is
#' an error. When only a few individuals are of interest, subset before
#' merging - both to keep the result small and to sidestep duplicate IDs
#' elsewhere in a large file.
#'
#' When `loci` is left at its default and the two inputs' panels differ, a
#' message reports how many loci are unique to each side and whether they
#' were kept or dropped - the two panels shown in the examples below are
#' typical, so this fires often and is meant to be seen once, not hunted for.
#'
#' @seealso [read_genotypes()], [as_geno_data()], [check_trio()].
#'
#' @examples
#' p_path <- system.file(
#'   "extdata", "Parents_FRH2024_broodstock_reference_genotypes.csv",
#'   package = "genoCompaR"
#' )
#' u_path <- system.file(
#'   "extdata", "Unknowns_CFS_WY26_salvage_v3_03262026.csv",
#'   package = "genoCompaR"
#' )
#'
#' parents <- read_genotypes(p_path, id_col = "SAMPLE_ID", group_col = NULL)
#' # the unknowns file contains a few repeated IDs, so keep the row of
#' # interest before building a geno_data object
#' raw <- utils::read.csv(u_path, colClasses = "character", check.names = FALSE)
#' offspring <- as_geno_data(
#'   raw[raw$SAMPLE_ID == "C260383CVP", ],
#'   id_col = "SAMPLE_ID", group_col = NULL
#' )
#'
#' trio <- merge_genotypes(
#'   offspring,
#'   parents[parents$SAMPLE_ID %in% c("C240716FRH", "C240722FRH"), ]
#' )
#' check_trio(trio, "C260383CVP", "C240716FRH", "C240722FRH")
#'
#' @export
merge_genotypes <- function(x, y, id_col = NULL, loci = NULL,
                             keep_unmatched_loci = FALSE) {
  if (!inherits(x, "geno_data") || !inherits(y, "geno_data")) {
    stop(
      "`x` and `y` must both be geno_data objects (see read_genotypes() ",
      "or as_geno_data()).",
      call. = FALSE
    )
  }

  shared <- intersect(get_loci(x), get_loci(y))
  candidate_pool <- if (keep_unmatched_loci) {
    union(get_loci(x), get_loci(y))
  } else {
    shared
  }

  if (is.null(loci)) {
    loci <- candidate_pool
    x_only <- setdiff(get_loci(x), get_loci(y))
    y_only <- setdiff(get_loci(y), get_loci(x))
    if (length(x_only) > 0 || length(y_only) > 0) {
      disposition <- if (keep_unmatched_loci) {
        "kept, NA-filled for the input that lacks them (keep_unmatched_loci = TRUE)"
      } else {
        "dropped as not shared by both inputs (keep_unmatched_loci = TRUE keeps them)"
      }
      message(sprintf(
        "merge_genotypes(): %d locus/loci only in `x` and %d only in `y` %s.",
        length(x_only), length(y_only), disposition
      ))
    }
  } else {
    unknown <- setdiff(loci, candidate_pool)
    if (length(unknown) > 0) {
      stop(
        "Locus/loci not present in ",
        if (keep_unmatched_loci) "either input: " else "both inputs: ",
        paste(utils::head(unknown, 5), collapse = ", "),
        if (length(unknown) > 5) ", ..." else "",
        call. = FALSE
      )
    }
  }
  if (length(loci) == 0) {
    stop("`x` and `y` share no loci.", call. = FALSE)
  }

  if (is.null(id_col)) id_col <- .id_col(x)
  locus_cols <- as.vector(rbind(paste0(loci, ".1"), paste0(loci, ".2")))

  x_gc <- .group_col(x)
  y_gc <- .group_col(y)
  keep_group <- !is.null(x_gc) && !is.null(y_gc)

  take <- function(d, gc) {
    df <- as.data.frame(d)
    absent <- setdiff(locus_cols, colnames(df))
    if (length(absent) > 0) df[absent] <- NA_character_

    out <- df[, c(.id_col(d), if (keep_group) gc, locus_cols), drop = FALSE]
    names(out)[1] <- id_col
    if (keep_group) names(out)[2] <- "group"
    out
  }
  combined <- rbind(take(x, x_gc), take(y, y_gc))

  as_geno_data(
    combined,
    id_col = id_col,
    group_col = if (keep_group) "group" else NULL,
    loci = loci
  )
}

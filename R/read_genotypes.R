#' Read a wide-format genotype/haplotype file
#'
#' Reads a CSV genotype/haplotype call sheet in "wide" format, where each
#' locus is represented by two adjacent allele columns named `<locus>.1` and
#' `<locus>.2` (as exported by many GT-seq / haplotype calling pipelines),
#' and returns a [tibble::tibble] of class `geno_data` carrying attributes
#' that the rest of genoCompaR relies on: the sample ID column, an optional
#' group column, and the list of locus names.
#'
#' @param file Path to a CSV file, or any input accepted by
#'   [utils::read.csv()].
#' @param id_col Name of the column holding unique sample identifiers.
#'   Defaults to `"indiv.ID"`.
#' @param group_col Name of an optional grouping column (e.g. population,
#'   run, or study group). Defaults to `"group"`; set to `NULL` (or the name
#'   of a column not present in the file) if there is no such column.
#' @param na_strings Character vector of strings treated as missing allele
#'   calls. Defaults to `c("NA", "")`.
#'
#' @return A `geno_data` object: a [tibble::tibble] with one row per sample,
#'   the `id_col` and (if present) `group_col` columns, and two columns per
#'   locus (`<locus>.1`, `<locus>.2`). Carries attributes `"loci"` (character
#'   vector of locus names), `"id_col"`, and `"group_col"`.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' geno
#' get_loci(geno)[1:5]
#'
#' @export
read_genotypes <- function(file,
                            id_col = "indiv.ID",
                            group_col = "group",
                            na_strings = c("NA", "")) {
  raw <- utils::read.csv(
    file,
    colClasses = "character",
    na.strings = na_strings,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!id_col %in% names(raw)) {
    stop(
      "id_col '", id_col, "' not found in ", file,
      ". Is your ID column named '", names(raw)[1], "'?",
      call. = FALSE
    )
  }
  if (anyDuplicated(raw[[id_col]])) {
    dupes <- unique(raw[[id_col]][duplicated(raw[[id_col]])])
    stop(
      "Duplicate sample IDs in id_col '", id_col, "': ",
      paste(utils::head(dupes, 5), collapse = ", "),
      if (length(dupes) > 5) ", ..." else "",
      call. = FALSE
    )
  }

  has_group <- !is.null(group_col) && group_col %in% names(raw)
  meta_cols <- c(id_col, if (has_group) group_col)

  locus_cols <- setdiff(names(raw), meta_cols)
  loci <- .validate_locus_columns(locus_cols)

  out <- tibble::as_tibble(raw[, c(meta_cols, locus_cols), drop = FALSE])
  structure(
    out,
    class = c("geno_data", class(out)),
    id_col = id_col,
    group_col = if (has_group) group_col else NA_character_,
    loci = loci
  )
}

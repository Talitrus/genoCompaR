#' Convert NMFS PBT/CKMR full-sibling family output to CKMRsim PO format
#'
#' Reads an NMFS parentage-based tagging (PBT) / CKMR full-sibling family
#' workbook (such as `Salvage2026_SanJoaquin_ckmr.xlsx`), where each row lists a
#' reconstructed family (`ma`, `pa`, `fs_fam_size`, and `kid1`, `kid2`, ...),
#' and expands every valid offspring-parent pair into the pairwise
#' parent-offspring (PO) table layout produced by CKMRsim (matching e.g.
#' `Coleman_BY25_PO_results.tsv`). Columns not present in the NMFS family
#' summary (`logl_ratio`, `num_loc`, `total_incompat`, `total_compat`,
#' `fraction_incompat`) are populated with `NA`.
#'
#' Summary or total rows at the bottom of the sheet (where `ma` and `pa` or all
#' `kid*` columns are empty / `NA`) are automatically ignored.
#'
#' @param file Path to an Excel workbook (`.xlsx` or `.xls`), or a
#'   [data.frame] / [tibble::tibble] already containing `ma`, `pa`, and one or
#'   more `kid*` columns.
#' @param out_file Optional path to write the output as a tab-separated (`.tsv`)
#'   file. When `NULL` (the default), no file is written and the tibble is
#'   returned directly.
#' @param sheet Sheet to read when `file` is an Excel workbook path: either a
#'   sheet name (character) or 1-based sheet index (numeric). Defaults to `1`.
#' @param strip_kid_suffix Logical (default `TRUE`). When `TRUE`, removes a
#'   trailing `_<ROSA/haplotype>` suffix from offspring (`D2_indiv`) sample IDs
#'   (e.g. `"C260698CVP_LNNLLNLLLLLL"` becomes `"C260698CVP"`), while preserving
#'   parent (`D1_indiv`) IDs verbatim.
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{`D2_indiv`}{Character. Offspring (juvenile / kid) sample ID.}
#'     \item{`D1_indiv`}{Character. Candidate parent (`ma` or `pa`) sample ID.}
#'     \item{`logl_ratio`}{Numeric (`NA_real_`). Log-likelihood ratio placeholder.}
#'     \item{`num_loc`}{Integer (`NA_integer_`). Number of loci compared placeholder.}
#'     \item{`total_incompat`}{Integer (`NA_integer_`). Incompatible loci count placeholder.}
#'     \item{`total_compat`}{Integer (`NA_integer_`). Compatible loci count placeholder.}
#'     \item{`fraction_incompat`}{Numeric (`NA_real_`). Fraction incompatible placeholder.}
#'   }
#'
#' @examples
#' df <- data.frame(
#'   index = c(1, 2, NA),
#'   ma = c("C250002MER", "T166266_7949_scarfF_2018", NA),
#'   pa = c("NA", "T246584_SJRRP2025ADULT132_returnM", NA),
#'   fs_fam_size = c(2, 1, 3),
#'   kid1 = c("C260698CVP_LNNLLNLLLLLL", "C259512CVP_ENNEENEEEEEE", NA),
#'   kid2 = c("C260734CVP_LNNLLNLLLLLL", NA, NA),
#'   stringsAsFactors = FALSE
#' )
#' nmfs_pbt2ckmrout(df)
#'
#' @export
nmfs_pbt2ckmrout <- function(file,
                             out_file = NULL,
                             sheet = 1,
                             strip_kid_suffix = TRUE) {
  if (!is.logical(strip_kid_suffix) || length(strip_kid_suffix) != 1L || is.na(strip_kid_suffix)) {
    stop("`strip_kid_suffix` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(out_file)) {
    if (!is.character(out_file) || length(out_file) != 1L || is.na(out_file) || !nzchar(out_file)) {
      stop("`out_file` must be a single non-empty file path or NULL.", call. = FALSE)
    }
  }

  if (is.data.frame(file)) {
    raw <- file
  } else if (is.character(file) && length(file) == 1L && !is.na(file) && nzchar(file)) {
    if (!file.exists(file)) {
      stop("File not found: ", file, call. = FALSE)
    }
    raw <- readxl::read_excel(
      path = file,
      sheet = sheet,
      col_types = "text",
      na = c("", "NA"),
      .name_repair = "minimal"
    )
  } else {
    stop("`file` must be a path to an Excel file or a data.frame.", call. = FALSE)
  }

  col_names <- names(raw)
  ma_idx <- which(tolower(col_names) == "ma")
  pa_idx <- which(tolower(col_names) == "pa")
  if (length(ma_idx) != 1L || length(pa_idx) != 1L) {
    stop("Input must contain 'ma' and 'pa' columns.", call. = FALSE)
  }

  kid_idx <- grep("^kid\\d*$", col_names, ignore.case = TRUE)
  if (length(kid_idx) == 0L) {
    stop(
      "Input must contain at least one 'kid' column (e.g. 'kid1', 'kid2', ...).",
      call. = FALSE
    )
  }

  d2_list <- list()
  d1_list <- list()
  pair_idx <- 0L

  for (i in seq_len(nrow(raw))) {
    ma_val <- trimws(as.character(raw[[ma_idx]][i]))
    pa_val <- trimws(as.character(raw[[pa_idx]][i]))

    parents <- c(
      if (.is_called(ma_val)) ma_val,
      if (.is_called(pa_val)) pa_val
    )
    if (length(parents) == 0L) next

    kids <- vapply(
      kid_idx,
      function(col_j) trimws(as.character(raw[[col_j]][i])),
      character(1),
      USE.NAMES = FALSE
    )
    kids <- kids[.is_called(kids)]
    if (length(kids) == 0L) next

    if (strip_kid_suffix) {
      kids <- sub("_[^_]+$", "", kids)
    }

    for (k in kids) {
      for (p in parents) {
        pair_idx <- pair_idx + 1L
        d2_list[[pair_idx]] <- k
        d1_list[[pair_idx]] <- p
      }
    }
  }

  d2_out <- if (pair_idx > 0L) unlist(d2_list, use.names = FALSE) else character(0)
  d1_out <- if (pair_idx > 0L) unlist(d1_list, use.names = FALSE) else character(0)
  n <- length(d2_out)

  out <- tibble::tibble(
    D2_indiv = d2_out,
    D1_indiv = d1_out,
    logl_ratio = rep(NA_real_, n),
    num_loc = rep(NA_integer_, n),
    total_incompat = rep(NA_integer_, n),
    total_compat = rep(NA_integer_, n),
    fraction_incompat = rep(NA_real_, n)
  )

  if (!is.null(out_file)) {
    utils::write.table(
      out,
      file = out_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE,
      na = "NA"
    )
  }

  out
}


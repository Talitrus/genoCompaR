#' Compare replicate samples within a genotype dataset
#'
#' Compares multiple samples in a `geno_data` object that are putative
#' replicates (e.g. technical replicates, re-extracted or re-run samples,
#' or suspected duplicate individuals). Computes pairwise exact-match
#' genotype concordance and distance, performs a 2D ordination (PCoA)
#' to visually map relationships between samples, and identifies
#' "culprit" loci with discrepant calls across replicates.
#'
#' @details
#' Genotype distance is calculated as the proportional mismatch distance
#' across all compared loci (loci successfully called in both samples):
#' \deqn{d_{ij} = \frac{\text{mismatching loci}_{ij}}{\text{compared loci}_{ij}} = 1 - \frac{\text{matching loci}_{ij}}{\text{compared loci}_{ij}}}
#' At each locus, genotypes are evaluated as exact diploid matches
#' (both alleles matching, unordered). Self-distance \eqn{d_{ii}} is 0.
#'
#' Ordination is performed using Classical Multidimensional Scaling
#' (Principal Coordinates Analysis, PCoA) via [stats::cmdscale()] on the
#' pairwise distance matrix, projecting samples into 2D space where
#' Euclidean distance approximates genotype mismatch distance.
#'
#' Discrepant loci (where at least one replicate differs from the others)
#' are categorized into:
#' \describe{
#'   \item{`allelic_dropout`}{Calls share alleles with one or more samples
#'     being homozygous for an allele present in a heterozygote (e.g.
#'     `A/A` in one sample, `A/B` in another), typical of low depth or
#'     dropout.}
#'   \item{`discordant`}{Calls have conflicting alleles with no shared
#'     alleles (e.g. `A/A` vs `B/B` or `A/B` vs `C/D`), indicative of
#'     contamination or sample swap.}
#'   \item{`multiallelic_discordant`}{More than two distinct alleles observed
#'     across the replicates (e.g. `A/B` and `A/C`).}
#' }
#'
#' @param x A `geno_data` object, e.g. from [read_genotypes()].
#' @param samples Character vector of sample IDs (or numeric row indices)
#'   to compare. If `NULL` (the default), all samples in `x` are compared
#'   (assuming `x` has already been subsetted).
#' @param min_loci_compared Minimum number of loci both samples in a pair
#'   must have been genotyped at for that pair to be scored. Defaults to `10`.
#' @param loci Character vector of loci to use. Defaults to all loci in `x`.
#'
#' @return An object of class `geno_replicates` (a list) containing:
#' \describe{
#'   \item{`summary`}{A tibble of pairwise comparisons with columns
#'     `sample1`, `sample2`, `n_loci_compared`, `n_match`, `n_mismatch`,
#'     `pct_match`, and `distance`.}
#'   \item{`distance_matrix`}{A `dist` object and symmetric numeric matrix
#'     of pairwise mismatch distances.}
#'   \item{`discrepancies`}{A tibble of all loci where samples disagree,
#'     with locus name, number of typed samples, distinct genotype counts,
#'     discrepancy classification (`discrepancy_type`), and each sample's
#'     individual call.}
#'   \item{`ordination`}{A list with `$points` (2D coordinates) and
#'     `$var_explained` (percentage of variance explained by Axis 1 & 2).}
#'   \item{`samples`}{Character vector of sample IDs compared.}
#'   \item{`n_loci_total`}{Total number of loci evaluated.}
#' }
#'
#' @seealso [plot.geno_replicates()] to visualize sample ordination;
#'   [compare_genotypes()] for detailed two-sample comparisons;
#'   [find_duplicate_genotypes()] to scan an entire panel for matches.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#'
#' # Compare the first 4 samples
#' ids <- geno$indiv.ID[1:4]
#' reps <- compare_replicates(geno, samples = ids)
#' reps
#'
#' # Inspect pairwise concordance summary
#' reps$summary
#'
#' # View discrepant loci (if any)
#' reps$discrepancies
#'
#' # Plot ordination
#' plot(reps)
#'
#' @export
compare_replicates <- function(x,
                               samples = NULL,
                               min_loci_compared = 10,
                               loci = get_loci(x)) {
  if (!inherits(x, "geno_data")) {
    stop("`x` must be a geno_data object (see read_genotypes()).", call. = FALSE)
  }

  all_ids <- as.character(x[[.id_col(x)]])

  if (is.null(samples)) {
    sample_ids <- all_ids
    row_idx <- seq_len(nrow(x))
  } else if (is.numeric(samples)) {
    if (any(samples < 1 | samples > nrow(x))) {
      stop(
        "Numeric `samples` indices must be between 1 and ", nrow(x), ".",
        call. = FALSE
      )
    }
    sample_ids <- all_ids[samples]
    row_idx <- samples
  } else {
    sample_ids <- as.character(samples)
    missing_ids <- setdiff(sample_ids, all_ids)
    if (length(missing_ids) > 0) {
      stop(
        "Sample(s) not found in '", .id_col(x), "': ",
        paste(missing_ids, collapse = ", "),
        call. = FALSE
      )
    }
    if (anyDuplicated(sample_ids)) {
      warning("Duplicate sample IDs in `samples` were deduplicated.", call. = FALSE)
      sample_ids <- unique(sample_ids)
    }
    row_idx <- match(sample_ids, all_ids)
  }

  n <- length(sample_ids)
  if (n < 2) {
    stop("Need at least 2 samples to compare replicates.", call. = FALSE)
  }

  # Validate loci
  all_loci <- get_loci(x)
  bad_loci <- setdiff(loci, all_loci)
  if (length(bad_loci) > 0) {
    stop(
      "The following loci were not found in `x`: ",
      paste(utils::head(bad_loci, 5), collapse = ", "),
      if (length(bad_loci) > 5) ", ..." else "",
      call. = FALSE
    )
  }

  sub_x <- x[row_idx, , drop = FALSE]

  # 1. Pairwise counts & distance matrix
  counts <- .pairwise_counts(sub_x, loci)
  comp_cnt <- counts$compared_count
  match_cnt <- counts$match_count
  mismatch_cnt <- comp_cnt - match_cnt

  ok <- comp_cnt > 0 & comp_cnt >= min_loci_compared
  pct_match_mat <- ifelse(ok, 100 * match_cnt / comp_cnt, NA_real_)
  dist_mat <- ifelse(ok, mismatch_cnt / comp_cnt, NA_real_)
  diag(dist_mat) <- 0
  dimnames(dist_mat) <- list(sample_ids, sample_ids)
  dimnames(pct_match_mat) <- list(sample_ids, sample_ids)

  dist_obj <- stats::as.dist(dist_mat)

  # 2. Pairwise summary tibble
  pair_idx <- which(upper.tri(dist_mat), arr.ind = TRUE)
  order_idx <- order(pair_idx[, 1], pair_idx[, 2])
  pair_idx <- pair_idx[order_idx, , drop = FALSE]

  summary_tbl <- tibble::tibble(
    sample1 = sample_ids[pair_idx[, 1]],
    sample2 = sample_ids[pair_idx[, 2]],
    n_loci_compared = comp_cnt[pair_idx],
    n_match = match_cnt[pair_idx],
    n_mismatch = mismatch_cnt[pair_idx],
    pct_match = pct_match_mat[pair_idx],
    distance = dist_mat[pair_idx]
  )

  gc <- .group_col(x)
  if (!is.null(gc)) {
    summary_tbl$group1 <- sub_x[[gc]][pair_idx[, 1]]
    summary_tbl$group2 <- sub_x[[gc]][pair_idx[, 2]]
  }

  # 3. Discrepancy & culprit locus analysis
  disc_loci <- character(0)
  disc_ntyped <- integer(0)
  disc_ndistinct <- integer(0)
  disc_genotypes <- character(0)
  disc_types <- character(0)
  disc_calls_list <- vector("list", n)
  for (s in seq_len(n)) {
    disc_calls_list[[s]] <- character(0)
  }

  for (li in seq_along(loci)) {
    locus <- loci[li]
    c1 <- sub_x[[paste0(locus, ".1")]]
    c2 <- sub_x[[paste0(locus, ".2")]]
    keys <- .canon_genotype_vec(c1, c2)

    valid_idx <- which(!is.na(keys))
    valid_keys <- keys[valid_idx]
    uniq_calls <- unique(valid_keys)

    if (length(uniq_calls) > 1) {
      disc_loci <- c(disc_loci, locus)
      disc_ntyped <- c(disc_ntyped, length(valid_keys))
      disc_ndistinct <- c(disc_ndistinct, length(uniq_calls))

      # Classify discrepancy
      allele_list <- lapply(strsplit(uniq_calls, "/"), unique)
      all_alleles <- unique(unlist(allele_list))

      if (length(all_alleles) <= 2) {
        # Biallelic: check if heterozygous call is present alongside homozygote(s)
        has_het <- any(vapply(allele_list, function(a) length(a) == 2, logical(1)))
        discrep_type <- if (has_het) "allelic_dropout" else "discordant"
      } else {
        # More than two distinct alleles across replicates
        discrep_type <- "multiallelic_discordant"
      }
      disc_types <- c(disc_types, discrep_type)

      tab <- table(valid_keys)
      tab_str <- paste(sprintf("%s (%d)", names(tab), as.integer(tab)), collapse = ", ")
      disc_genotypes <- c(disc_genotypes, tab_str)

      for (s in seq_len(n)) {
        call_val <- if (is.na(keys[s])) "NA" else keys[s]
        disc_calls_list[[s]] <- c(disc_calls_list[[s]], call_val)
      }
    }
  }

  discrep_df <- tibble::tibble(
    locus = disc_loci,
    n_typed = disc_ntyped,
    n_distinct_genotypes = disc_ndistinct,
    distinct_genotypes = disc_genotypes,
    discrepancy_type = disc_types
  )
  for (s in seq_len(n)) {
    discrep_df[[paste0("call_", sample_ids[s])]] <- disc_calls_list[[s]]
  }

  if (nrow(discrep_df) > 0) {
    discrep_df <- discrep_df[order(-discrep_df$n_distinct_genotypes, discrep_df$locus), ]
  }

  # 4. Ordination (PCoA via cmdscale)
  points <- NULL
  var_exp <- c(NA_real_, NA_real_)

  if (n == 2) {
    # 2 points collapse to 1D
    d_val <- dist_mat[1, 2]
    if (!is.na(d_val)) {
      points <- matrix(
        c(-d_val / 2, d_val / 2, 0, 0),
        nrow = 2, ncol = 2,
        dimnames = list(sample_ids, c("PCoA1", "PCoA2"))
      )
      var_exp <- c(100, 0)
    }
  } else if (n >= 3) {
    if (any(is.na(dist_mat[upper.tri(dist_mat)]))) {
      warning(
        "Pairwise distance matrix contains NA values due to pairs below `min_loci_compared`. ",
        "Ordination could not be computed.",
        call. = FALSE
      )
    } else {
      fit <- stats::cmdscale(dist_obj, k = 2, eig = TRUE)
      points <- fit$points
      colnames(points) <- c("PCoA1", "PCoA2")
      rownames(points) <- sample_ids

      pos_eig <- fit$eig[fit$eig > 0]
      total_var <- sum(pos_eig)
      if (total_var > 0 && length(fit$eig) >= 2) {
        var_exp <- c(max(0, fit$eig[1]), max(0, fit$eig[2])) / total_var * 100
      }
    }
  }

  structure(
    list(
      summary = summary_tbl,
      distance_matrix = dist_obj,
      discrepancies = discrep_df,
      ordination = list(
        points = points,
        var_explained = var_exp
      ),
      samples = sample_ids,
      n_loci_total = length(loci)
    ),
    class = "geno_replicates"
  )
}

#' Plot replicate sample ordination
#'
#' Visualizes the 2D PCoA ordination of replicate samples from a
#' `geno_replicates` object returned by [compare_replicates()].
#'
#' @param x A `geno_replicates` object from [compare_replicates()].
#' @param engine One of `"base"` (default) or `"ggplot2"`. If `"ggplot2"` is
#'   selected, the package `ggplot2` must be installed.
#' @param labels Logical; if `TRUE` (default), prints sample ID labels near points.
#' @param ... Additional arguments passed to base [graphics::plot()] or ignored.
#'
#' @return Invisibly returns the 2D coordinate matrix (for `engine = "base"`)
#'   or the `ggplot` object (for `engine = "ggplot2"`).
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' reps <- compare_replicates(geno, samples = geno$indiv.ID[1:4])
#'
#' # Base graphics (default)
#' plot(reps)
#'
#' \dontrun{
#' # Optional ggplot2 output
#' plot(reps, engine = "ggplot2")
#' }
#'
#' @export
plot.geno_replicates <- function(x,
                                 engine = c("base", "ggplot2"),
                                 labels = TRUE,
                                 ...) {
  engine <- match.arg(engine)

  pts <- x$ordination$points
  if (is.null(pts)) {
    stop("No ordination coordinates available to plot.", call. = FALSE)
  }

  ve <- x$ordination$var_explained
  xlab <- if (!is.na(ve[1])) sprintf("PCoA 1 (%.1f%%)", ve[1]) else "PCoA 1"
  ylab <- if (!is.na(ve[2])) sprintf("PCoA 2 (%.1f%%)", ve[2]) else "PCoA 2"
  title_str <- sprintf("Replicate Ordination (n = %d)", nrow(pts))

  if (engine == "base") {
    xr <- range(pts[, 1])
    yr <- range(pts[, 2])
    xd <- diff(xr)
    if (xd == 0) xd <- 0.1
    yd <- diff(yr)
    if (yd == 0) yd <- 0.1

    xlim <- c(xr[1] - 0.20 * xd, xr[2] + 0.20 * xd)
    ylim <- c(yr[1] - 0.20 * yd, yr[2] + 0.20 * yd)

    graphics::plot(
      pts[, 1], pts[, 2],
      xlab = xlab, ylab = ylab,
      main = title_str,
      xlim = xlim, ylim = ylim,
      pch = 19, col = "#2C3E50", cex = 1.3,
      ...
    )
    graphics::abline(h = 0, v = 0, lty = 2, col = "grey80")

    if (labels) {
      graphics::text(
        pts[, 1], pts[, 2],
        labels = rownames(pts),
        pos = 3, cex = 0.85, col = "#1A252F"
      )
    }

    invisible(pts)
  } else {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop(
        "Package 'ggplot2' is required for `engine = 'ggplot2'`. ",
        "Please install it with install.packages('ggplot2') or use `engine = 'base'`.",
        call. = FALSE
      )
    }

    df <- data.frame(
      sample = rownames(pts),
      PCoA1 = pts[, 1],
      PCoA2 = pts[, 2],
      stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PCoA1, y = .data$PCoA2)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
      ggplot2::geom_point(color = "#2C3E50", size = 3) +
      ggplot2::labs(
        title = title_str,
        x = xlab,
        y = ylab
      ) +
      ggplot2::theme_minimal()

    if (labels) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$sample),
        vjust = -1, size = 3.5
      )
    }

    if (interactive()) {
      print(p)
    }
    invisible(p)
  }
}

#' @rdname plot.geno_replicates
#' @export
plot_replicates <- function(x,
                            engine = c("base", "ggplot2"),
                            labels = TRUE,
                            ...) {
  plot.geno_replicates(x, engine = engine, labels = labels, ...)
}

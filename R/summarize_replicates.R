#' Summarize trends across multiple replicate comparisons
#'
#' Takes multiple outputs from [compare_replicates()] (e.g. replicate sets
#' from different individuals or runs) and summarizes them to detect
#' systematic panel-wide trends: identifying "chronic culprit" loci that fail
#' repeatedly across different individuals, and flagging outlier replicate
#' sets that fall below expected concordance (potential sample swaps or
#' mix-ups).
#'
#' @details
#' In technical replicate or duplicate genotyping, true replicates typically
#' exhibit high concordance (>= 95%). This function aggregates across multiple
#' replicate comparison runs to surface two key QC signals:
#'
#' 1. **Replicate Set Outliers:** Any replicate set whose minimum or mean
#'    concordance dips below `min_concordance` (default 95%) is flagged as a
#'    potential sample mix-up, mislabeling, or cross-contamination.
#' 2. **Chronic Culprit Loci:** Loci that show discrepancies across multiple
#'    independent replicate sets (exceeding `locus_error_threshold`) are
#'    flagged as systematic error culprits (e.g. primer site mutations,
#'    consistently poor amplification, or paralogous interference).
#'
#' The function records observed discrepancy classifications (`allelic_dropout`,
#' `discordant`, `multiallelic_discordant`) for each locus to assist in
#' diagnosing whether differences look like allele drop-outs (sharing an allele)
#' or conflicting alleles.
#'
#' @param replicates A list of `geno_replicates` objects from
#'   [compare_replicates()], or a single `geno_replicates` object.
#' @param ... Additional `geno_replicates` objects if passed individually
#'   rather than as a list.
#' @param min_concordance Minimum percentage match (0-100) expected for
#'   replicate pairs. Sets with concordance below this threshold are flagged
#'   as outliers (`flagged_outlier = TRUE`). Defaults to `95`.
#' @param locus_error_threshold Minimum proportion (0-1) of evaluated replicate
#'   sets in which a locus can be discrepant before being flagged as a chronic
#'   culprit. Defaults to `0.10` (10% of sets).
#'
#' @return An object of class `replicate_trends` (a list) containing:
#' \describe{
#'   \item{`locus_trends`}{A tibble summarizing every evaluated locus across
#'     all replicate sets, giving `n_sets_evaluated`, `n_sets_discrepant`,
#'     `discrepancy_rate` (percentage of sets), counts of `n_dropout`,
#'     `n_discordant`, and `n_multiallelic`, and a `flagged` boolean.}
#'   \item{`set_summary`}{A tibble summarizing each replicate set:
#'     `set_name`, `n_samples`, `n_pairs`, `min_pct_match`, `mean_pct_match`,
#'     `median_pct_match`, `mean_distance`, `n_discrepant_loci`, and
#'     `flagged_outlier`.}
#'   \item{`global`}{A list of panel-wide statistics: total sets, total samples,
#'     total pairwise comparisons, mean/median concordance, clean locus rate,
#'     and count of flagged loci and sets.}
#'   \item{`all_discrepancies`}{A tidy tibble combining all discrepant locus
#'     records across all replicate sets, with a `replicate_set` column.}
#'   \item{`config`}{List of parameter thresholds used.}
#' }
#'
#' @seealso [compare_replicates()] to generate replicate comparison objects;
#'   [plot.replicate_trends()] to visualize chronic culprits or concordance.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#'
#' # Compare two different sets of replicate samples
#' set1 <- compare_replicates(geno, samples = geno$indiv.ID[1:3])
#' set2 <- compare_replicates(geno, samples = geno$indiv.ID[4:6])
#'
#' # Summarize trends across both sets
#' trends <- summarize_replicates(list(Fish1 = set1, Fish2 = set2))
#' trends
#'
#' # Inspect chronic culprit loci
#' trends$locus_trends
#'
#' # Inspect replicate set summaries
#' trends$set_summary
#'
#' # Plot top culprit loci
#' plot(trends, type = "culprits")
#'
#' @export
summarize_replicates <- function(replicates,
                                 ...,
                                 min_concordance = 95,
                                 locus_error_threshold = 0.10) {
  if (missing(replicates)) {
    reps_list <- list(...)
  } else if (inherits(replicates, "geno_replicates")) {
    reps_list <- list(replicates, ...)
  } else if (is.list(replicates)) {
    reps_list <- c(replicates, list(...))
  } else {
    stop("`replicates` must be a list of `geno_replicates` objects (see compare_replicates()).", call. = FALSE)
  }

  if (length(reps_list) == 0) {
    stop("Need at least one `geno_replicates` object to summarize.", call. = FALSE)
  }

  is_gr <- vapply(reps_list, function(obj) inherits(obj, "geno_replicates"), logical(1))
  if (!all(is_gr)) {
    stop("All inputs must be `geno_replicates` objects (from compare_replicates()).", call. = FALSE)
  }

  if (!is.numeric(min_concordance) || min_concordance < 0 || min_concordance > 100) {
    stop("`min_concordance` must be a percentage between 0 and 100.", call. = FALSE)
  }
  if (!is.numeric(locus_error_threshold) || locus_error_threshold < 0 || locus_error_threshold > 1) {
    stop("`locus_error_threshold` must be a proportion between 0 and 1.", call. = FALSE)
  }

  # Assign clean names to sets
  set_names <- names(reps_list)
  if (is.null(set_names)) {
    set_names <- paste0("set_", seq_along(reps_list))
  } else {
    empty_idx <- which(set_names == "" | is.na(set_names))
    if (length(empty_idx) > 0) {
      set_names[empty_idx] <- paste0("set_", empty_idx)
    }
  }
  names(reps_list) <- set_names

  n_sets <- length(reps_list)

  # 1. Per-set summary
  set_rows <- vector("list", n_sets)
  for (i in seq_len(n_sets)) {
    rep_obj <- reps_list[[i]]
    pcts <- rep_obj$summary$pct_match[!is.na(rep_obj$summary$pct_match)]
    dists <- rep_obj$summary$distance[!is.na(rep_obj$summary$distance)]

    min_pct <- if (length(pcts) > 0) min(pcts) else NA_real_
    mean_pct <- if (length(pcts) > 0) mean(pcts) else NA_real_
    med_pct <- if (length(pcts) > 0) stats::median(pcts) else NA_real_
    mean_dist <- if (length(dists) > 0) mean(dists) else NA_real_

    n_disc <- nrow(rep_obj$discrepancies)
    flagged <- !is.na(min_pct) && (min_pct < min_concordance)

    set_rows[[i]] <- list(
      set_name = set_names[i],
      n_samples = length(rep_obj$samples),
      n_pairs = nrow(rep_obj$summary),
      min_pct_match = min_pct,
      mean_pct_match = mean_pct,
      median_pct_match = med_pct,
      mean_distance = mean_dist,
      n_discrepant_loci = n_disc,
      flagged_outlier = flagged
    )
  }
  set_summary <- tibble::as_tibble(do.call(rbind.data.frame, c(set_rows, stringsAsFactors = FALSE)))

  # 2. Combined discrepancies
  disc_dfs <- vector("list", n_sets)
  for (i in seq_len(n_sets)) {
    df <- reps_list[[i]]$discrepancies
    if (nrow(df) > 0) {
      base_cols <- c("locus", "n_typed", "n_distinct_genotypes", "distinct_genotypes", "discrepancy_type")
      sub_df <- df[, base_cols, drop = FALSE]
      sub_df$replicate_set <- set_names[i]
      disc_dfs[[i]] <- sub_df[, c("replicate_set", base_cols)]
    }
  }
  all_disc_raw <- do.call(rbind, disc_dfs)
  if (is.null(all_disc_raw) || nrow(all_disc_raw) == 0) {
    all_discrepancies <- tibble::tibble(
      replicate_set = character(0),
      locus = character(0),
      n_typed = integer(0),
      n_distinct_genotypes = integer(0),
      distinct_genotypes = character(0),
      discrepancy_type = character(0)
    )
  } else {
    all_discrepancies <- tibble::as_tibble(all_disc_raw)
  }

  # 3. Locus trends
  # Determine all evaluated loci
  all_loci_list <- lapply(reps_list, function(r) {
    if (!is.null(r$loci)) r$loci else unique(r$discrepancies$locus)
  })
  unique_loci <- unique(unlist(all_loci_list))

  if (length(unique_loci) == 0 && nrow(all_discrepancies) > 0) {
    unique_loci <- unique(all_discrepancies$locus)
  }

  L <- length(unique_loci)
  if (L > 0) {
    n_sets_eval <- integer(L)
    n_sets_disc <- integer(L)
    n_dropout <- integer(L)
    n_discordant <- integer(L)
    n_multiallelic <- integer(L)

    for (li in seq_len(L)) {
      loc <- unique_loci[li]
      # Count sets that evaluated this locus
      n_sets_eval[li] <- sum(vapply(all_loci_list, function(locs) loc %in% locs, logical(1)))
      # Discrepancies for this locus
      loc_disc <- all_discrepancies[all_discrepancies$locus == loc, ]
      n_sets_disc[li] <- nrow(loc_disc)
      n_dropout[li] <- sum(loc_disc$discrepancy_type == "allelic_dropout")
      n_discordant[li] <- sum(loc_disc$discrepancy_type == "discordant")
      n_multiallelic[li] <- sum(loc_disc$discrepancy_type == "multiallelic_discordant")
    }

    disc_rate <- ifelse(n_sets_eval > 0, 100 * n_sets_disc / n_sets_eval, 0)
    flagged_locus <- (disc_rate / 100) >= locus_error_threshold

    locus_trends <- tibble::tibble(
      locus = unique_loci,
      n_sets_evaluated = n_sets_eval,
      n_sets_discrepant = n_sets_disc,
      discrepancy_rate = disc_rate,
      n_dropout = n_dropout,
      n_discordant = n_discordant,
      n_multiallelic = n_multiallelic,
      flagged = flagged_locus
    )
    # Sort by discrepancy_rate desc, then n_sets_discrepant desc, then locus
    locus_trends <- locus_trends[order(-locus_trends$discrepancy_rate, -locus_trends$n_sets_discrepant, locus_trends$locus), ]
  } else {
    locus_trends <- tibble::tibble(
      locus = character(0),
      n_sets_evaluated = integer(0),
      n_sets_discrepant = integer(0),
      discrepancy_rate = numeric(0),
      n_dropout = integer(0),
      n_discordant = integer(0),
      n_multiallelic = integer(0),
      flagged = logical(0)
    )
  }

  # 4. Global panel metrics
  all_pairs_list <- lapply(seq_len(n_sets), function(i) {
    df <- reps_list[[i]]$summary
    df$replicate_set <- set_names[i]
    df
  })
  all_pairs <- tibble::as_tibble(do.call(rbind, all_pairs_list))

  all_pcts <- all_pairs$pct_match[!is.na(all_pairs$pct_match)]
  n_clean_loci <- sum(locus_trends$n_sets_discrepant == 0)
  pct_clean <- if (L > 0) 100 * n_clean_loci / L else 100

  global_metrics <- list(
    n_sets = n_sets,
    n_samples_total = sum(set_summary$n_samples),
    n_pairs_total = nrow(all_pairs),
    mean_concordance = if (length(all_pcts) > 0) mean(all_pcts) else NA_real_,
    median_concordance = if (length(all_pcts) > 0) stats::median(all_pcts) else NA_real_,
    n_loci_total = L,
    n_loci_clean = n_clean_loci,
    pct_loci_clean = pct_clean,
    n_sets_flagged = sum(set_summary$flagged_outlier),
    n_loci_flagged = sum(locus_trends$flagged)
  )

  structure(
    list(
      locus_trends = locus_trends,
      set_summary = set_summary,
      global = global_metrics,
      all_discrepancies = all_discrepancies,
      config = list(
        min_concordance = min_concordance,
        locus_error_threshold = locus_error_threshold
      )
    ),
    class = "replicate_trends"
  )
}

#' Plot replicate trends
#'
#' Visualizes either chronic culprit loci or concordance distributions from a
#' `replicate_trends` object created by [summarize_replicates()].
#'
#' @param x A `replicate_trends` object from [summarize_replicates()].
#' @param type One of `"culprits"` (default) to show loci with the highest
#'   discrepancy rates across replicate sets, or `"concordance"` to show the
#'   distribution of concordance rates across replicate sets.
#' @param engine One of `"base"` (default) or `"ggplot2"`. If `"ggplot2"` is
#'   selected, the package `ggplot2` must be installed.
#' @param top_n Maximum number of top culprit loci to show when
#'   `type = "culprits"`. Defaults to `15`.
#' @param ... Additional arguments passed to base graphics functions or ignored.
#'
#' @return Invisibly returns the plotted data table (or `ggplot` object if
#'   `engine = "ggplot2"`).
#'
#' @examples
#' path <- system.file(
#'   "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' geno <- read_genotypes(path)
#' set1 <- compare_replicates(geno, samples = geno$indiv.ID[1:3])
#' set2 <- compare_replicates(geno, samples = geno$indiv.ID[4:6])
#' trends <- summarize_replicates(list(Fish1 = set1, Fish2 = set2))
#'
#' # Plot top culprit loci with base graphics
#' plot(trends, type = "culprits")
#'
#' # Plot concordance distribution
#' plot(trends, type = "concordance")
#'
#' \dontrun{
#' # With ggplot2
#' plot(trends, type = "culprits", engine = "ggplot2")
#' }
#'
#' @export
plot.replicate_trends <- function(x,
                                  type = c("culprits", "concordance"),
                                  engine = c("base", "ggplot2"),
                                  top_n = 15,
                                  ...) {
  type <- match.arg(type)
  engine <- match.arg(engine)

  if (type == "culprits") {
    lt <- x$locus_trends
    disc_loci <- lt[lt$n_sets_discrepant > 0, , drop = FALSE]

    if (nrow(disc_loci) == 0) {
      message("No discrepant loci found across any replicate sets; nothing to plot.")
      return(invisible(disc_loci))
    }

    sub_lt <- utils::head(disc_loci, top_n)

    if (engine == "base") {
      # Horizontal barplot of discrepancy rate
      vals <- rev(sub_lt$discrepancy_rate)
      names_arg <- rev(sub_lt$locus)

      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(oldpar))
      graphics::par(mar = c(5, 12, 4, 2) + 0.1)

      max_val <- max(vals, na.rm = TRUE)
      xlim_max <- min(100, max_val * 1.25 + 5)

      bp <- graphics::barplot(
        vals,
        horiz = TRUE,
        names.arg = names_arg,
        las = 1,
        xlab = "% Replicate Sets Discrepant",
        main = sprintf("Top Culprit Loci across Replicate Sets (n = %d)", x$global$n_sets),
        col = "#2C3E50",
        border = NA,
        xlim = c(0, xlim_max),
        ...
      )
      graphics::grid(nx = NULL, ny = NA, col = "grey85", lty = "dotted")
      graphics::barplot(
        vals,
        horiz = TRUE,
        names.arg = names_arg,
        las = 1,
        col = "#2C3E50",
        border = NA,
        add = TRUE
      )
      graphics::text(
        x = vals,
        y = bp,
        labels = sprintf(" %.1f%%", vals),
        pos = 4,
        cex = 0.8
      )
      invisible(sub_lt)
    } else {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for `engine = 'ggplot2'`. Please install it or use `engine = 'base'`.", call. = FALSE)
      }

      sub_lt$locus_factor <- stats::reorder(factor(sub_lt$locus), sub_lt$discrepancy_rate)

      p <- ggplot2::ggplot(sub_lt, ggplot2::aes(x = .data$locus_factor, y = .data$discrepancy_rate)) +
        ggplot2::geom_col(fill = "#2C3E50", width = 0.7) +
        ggplot2::coord_flip() +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.1f%%", .data$discrepancy_rate)),
          hjust = -0.15, size = 3.2
        ) +
        ggplot2::labs(
          title = sprintf("Top Culprit Loci across Replicate Sets (n = %d)", x$global$n_sets),
          x = "Locus",
          y = "% Replicate Sets Discrepant"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::expand_limits(y = max(sub_lt$discrepancy_rate) * 1.2)

      if (interactive()) {
        print(p)
      }
      invisible(p)
    }
  } else {
    # type == "concordance"
    ss <- x$set_summary
    thresh <- x$config$min_concordance

    if (engine == "base") {
      vals <- ss$min_pct_match
      names_arg <- ss$set_name
      cols <- ifelse(ss$flagged_outlier, "#E74C3C", "#2980B9")

      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(oldpar))

      ylim_min <- max(0, min(vals, thresh, na.rm = TRUE) - 5)
      graphics::plot(
        seq_along(vals), vals,
        pch = 19, col = cols, cex = 1.3,
        ylim = c(ylim_min, 101),
        xlab = "Replicate Set",
        ylab = "Minimum Concordance (% match)",
        main = "Concordance across Replicate Sets",
        xaxt = "n",
        ...
      )
      graphics::axis(1, at = seq_along(vals), labels = names_arg, las = 2, cex.axis = 0.8)
      graphics::abline(h = thresh, lty = 2, col = "#E74C3C", lwd = 1.5)
      graphics::text(
        x = 1, y = thresh,
        labels = sprintf("Threshold (%g%%)", thresh),
        pos = 3, col = "#E74C3C", cex = 0.8
      )
      invisible(ss)
    } else {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for `engine = 'ggplot2'`. Please install it or use `engine = 'base'`.", call. = FALSE)
      }

      p <- ggplot2::ggplot(ss, ggplot2::aes(x = .data$set_name, y = .data$min_pct_match, color = .data$flagged_outlier)) +
        ggplot2::geom_hline(yintercept = thresh, linetype = "dashed", color = "#E74C3C") +
        ggplot2::geom_point(size = 3) +
        ggplot2::scale_color_manual(
          values = c("FALSE" = "#2980B9", "TRUE" = "#E74C3C"),
          labels = c("FALSE" = "Concordant", "TRUE" = "Outlier (< threshold)"),
          name = "Status"
        ) +
        ggplot2::labs(
          title = "Concordance across Replicate Sets",
          x = "Replicate Set",
          y = "Minimum Concordance (% match)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

      if (interactive()) {
        print(p)
      }
      invisible(p)
    }
  }
}

#' @rdname plot.replicate_trends
#' @export
plot_replicate_trends <- function(x,
                                  type = c("culprits", "concordance"),
                                  engine = c("base", "ggplot2"),
                                  top_n = 15,
                                  ...) {
  plot.replicate_trends(x, type = type, engine = engine, top_n = top_n, ...)
}

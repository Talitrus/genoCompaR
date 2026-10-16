#' @export
print.geno_data <- function(x, ...) {
  cat(sprintf("<geno_data> %d samples x %d loci\n", nrow(x), length(get_loci(x))))
  gc <- .group_col(x)
  if (!is.null(gc)) {
    groups <- unique(x[[gc]])
    cat(sprintf(
      "Groups (%s): %s\n", gc,
      paste(utils::head(groups, 10), collapse = ", ")
    ))
  }
  NextMethod()
}

#' @export
print.geno_comparison <- function(x, ...) {
  # The report below describes exactly one pair, reading every field from row
  # 1. Several results row-bound together (the documented bulk workflow) must
  # therefore fall through to normal tibble printing rather than silently
  # rendering as the first pair only.
  if (nrow(x) != 1) {
    return(NextMethod())
  }
  detail <- locus_detail(x)
  mendelian <- isTRUE(attr(x, "mendelian"))
  rule <- strrep("-", 50)

  fmt_pct <- function(p) if (is.na(p)) "NA" else sprintf("%.1f%%", p)

  cat("Genotype comparison:", x$sample1[[1]], "vs", x$sample2[[1]], "\n")
  cat(rule, "\n", sep = "")
  cat(sprintf("Loci in panel:          %d\n", x$n_loci_total[[1]]))
  cat(sprintf("Loci compared:          %d (both samples genotyped)\n", x$n_loci_compared[[1]]))
  cat(sprintf("Loci unique to one:     %d (only one sample genotyped)\n", x$n_loci_unique[[1]]))
  cat(sprintf("Matching genotypes:     %d\n", x$n_match[[1]]))
  cat(sprintf("Mismatching genotypes:  %d\n", x$n_mismatch[[1]]))
  cat(sprintf("Percent match:          %s\n", fmt_pct(x$pct_match[[1]])))

  if (mendelian) {
    incompat_idx <- which(detail$mendelian_incompatible %in% TRUE)
    cat(rule, "\n", sep = "")
    cat(sprintf(
      "Mendelian incompatibilities: %d (%s of compared loci)\n",
      length(incompat_idx), fmt_pct(x$pct_mendelian_incompatible[[1]])
    ))
    if (length(incompat_idx) > 0) {
      show_n <- min(length(incompat_idx), 15)
      sub <- detail[incompat_idx[seq_len(show_n)], ]
      report <- data.frame(
        locus = sub$locus,
        sample1 = paste0(sub$sample1_a1, "/", sub$sample1_a2),
        sample2 = paste0(sub$sample2_a1, "/", sub$sample2_a2),
        stringsAsFactors = FALSE
      )
      print(report, row.names = FALSE)
      if (length(incompat_idx) > show_n) {
        cat("... and", length(incompat_idx) - show_n, "more\n")
      }
    }
  }

  cat(rule, "\n", sep = "")
  cat("Use locus_detail() on the returned object for the full per-locus table.\n")
  invisible(x)
}

#' @export
print.geno_duplicates <- function(x, ...) {
  cat(sprintf(
    "<geno_duplicates> %d flagged pair(s) among %d samples (%d loci, threshold = %s)\n",
    nrow(x), attr(x, "n_samples"), attr(x, "n_loci"),
    sprintf("%.0f%%", 100 * attr(x, "threshold"))
  ))
  NextMethod()
}

#' @export
print.geno_paralog_scan <- function(x, ...) {
  cat(sprintf(
    "<geno_paralog_scan> %d of %d loci scanned, %d flagged (excess > %.2f, min_n = %d)\n",
    nrow(x), attr(x, "n_loci"), sum(x$flagged), attr(x, "threshold"),
    as.integer(attr(x, "min_n"))
  ))
  skipped <- attr(x, "n_loci") - nrow(x)
  if (skipped > 0) {
    cat(sprintf(
      "%d locus/loci skipped: fewer than %d complete genotypes.\n",
      skipped, as.integer(attr(x, "min_n"))
    ))
  }
  NextMethod()
}

#' @export
print.geno_conformed <- function(x, ...) {
  report <- attr(x, "conform_report")
  cat(sprintf(
    "<geno_conformed> %d samples x %d columns | %d anomal%s\n",
    nrow(x), ncol(x), length(report$anomalies),
    if (length(report$anomalies) == 1) "y" else "ies"
  ))
  cat("Use conform_report() for the QC report, as_geno_data() to compare samples.\n")
  NextMethod()
}

#' @export
print.conform_report <- function(x, ...) {
  cat(x$text, sep = "\n")
  cat("\n")
  invisible(x)
}

#' @export
print.geno_haps <- function(x, ...) {
  report <- attr(x, "gen_haps_report")
  cat(sprintf(
    "<geno_haps> %d samples x %d loci | %d anomal%s\n",
    nrow(x), (ncol(x) - 2L) %/% 2L, length(report$anomalies),
    if (length(report$anomalies) == 1) "y" else "ies"
  ))
  cat("Use gen_haps_report() for the QC report, as_geno_data() to compare samples.\n")
  NextMethod()
}

#' @export
print.gen_haps_report <- function(x, ...) {
  cat(x$text, sep = "\n")
  cat("\n")
  invisible(x)
}

#' @export
print.geno_trio <- function(x, ...) {
  # see print.geno_comparison() - this report describes a single trio
  if (nrow(x) != 1) {
    return(NextMethod())
  }
  detail <- locus_detail(x)
  rule <- strrep("-", 50)
  fmt_pct <- function(p) if (is.na(p)) "NA" else sprintf("%.1f%%", p)

  cat(
    "Trio parentage check: offspring", x$offspring[[1]],
    "| parents", x$parent1[[1]], "&", x$parent2[[1]], "\n"
  )
  cat(rule, "\n", sep = "")
  cat(sprintf("Loci in panel:            %d\n", x$n_loci_total[[1]]))
  cat(sprintf("Loci compared:            %d (all three samples genotyped)\n", x$n_loci_compared[[1]]))
  cat(sprintf("Loci unique to one/two:   %d (at least one sample missing)\n", x$n_loci_unique[[1]]))
  cat(sprintf("Compatible loci:          %d\n", x$n_compatible[[1]]))
  cat(sprintf("Incompatible loci:        %d\n", x$n_incompatible[[1]]))
  cat(sprintf("Percent compatible:       %s\n", fmt_pct(x$pct_compatible[[1]])))

  incompat_idx <- which(detail$status == "incompatible")
  if (length(incompat_idx) > 0) {
    cat(rule, "\n", sep = "")
    show_n <- min(length(incompat_idx), 15)
    sub <- detail[incompat_idx[seq_len(show_n)], ]
    report <- data.frame(
      locus = sub$locus,
      offspring = paste0(sub$offspring_a1, "/", sub$offspring_a2),
      parent1 = paste0(sub$parent1_a1, "/", sub$parent1_a2),
      parent2 = paste0(sub$parent2_a1, "/", sub$parent2_a2),
      stringsAsFactors = FALSE
    )
    print(report, row.names = FALSE)
    if (length(incompat_idx) > show_n) {
      cat("... and", length(incompat_idx) - show_n, "more\n")
    }
  }

  cat(rule, "\n", sep = "")
  cat("Use locus_detail() on the returned object for the full per-locus table.\n")
  invisible(x)
}

# Arithmetic on a match matrix (`1 - mm / 100`, `mm / 2`, ...) produces
# something that is no longer a matrix of match percentages, so it must not
# keep the geno_match_matrix class or the now-stale n_compared attribute -
# otherwise the print header would report statistics that no longer describe
# the values.
#' @export
Ops.geno_match_matrix <- function(e1, e2) {
  strip <- function(v) {
    if (inherits(v, "geno_match_matrix")) {
      attr(v, "n_compared") <- NULL
      class(v) <- setdiff(class(v), "geno_match_matrix")
    }
    v
  }
  if (missing(e2)) {
    do.call(.Generic, list(strip(e1)))
  } else {
    do.call(.Generic, list(strip(e1), strip(e2)))
  }
}

#' @export
print.geno_match_matrix <- function(x, ...) {
  n <- nrow(x)
  cat(sprintf("<geno_match_matrix> %d x %d samples (values = %% loci matching)\n", n, n))

  off_diag <- x[upper.tri(x)]
  off_diag <- off_diag[!is.na(off_diag)]
  if (length(off_diag) > 0) {
    cat(sprintf(
      "Off-diagonal pct_match: min %.1f, median %.1f, max %.1f (%d pairs with data)\n",
      min(off_diag), stats::median(off_diag), max(off_diag), length(off_diag)
    ))
  }

  show_n <- min(n, 10)
  cat(sprintf("Showing first %d x %d samples (rounded):\n", show_n, show_n))
  print(round(unclass(x)[seq_len(show_n), seq_len(show_n), drop = FALSE], 1))
  if (n > show_n) {
    cat(sprintf(
      "... %d more sample(s) not shown. Index the matrix directly for the rest, e.g. `mm[1:20, 1:20]`.\n",
      n - show_n
    ))
  }
  invisible(x)
}

#' @export
print.geno_replicates <- function(x, ...) {
  n_samp <- length(x$samples)
  n_loci <- x$n_loci_total
  cat(sprintf("<geno_replicates> %d samples x %d loci\n", n_samp, n_loci))
  rule <- strrep("-", 50)
  cat(rule, "\n", sep = "")

  pcts <- x$summary$pct_match[!is.na(x$summary$pct_match)]
  dists <- x$summary$distance[!is.na(x$summary$distance)]

  if (length(pcts) > 0) {
    cat(sprintf(
      "Pairwise concordance (pct_match): min %.1f%%, median %.1f%%, max %.1f%%\n",
      min(pcts), stats::median(pcts), max(pcts)
    ))
    cat(sprintf(
      "Pairwise distance (mismatch rate): min %.3f, median %.3f, max %.3f\n",
      min(dists), stats::median(dists), max(dists)
    ))
  }

  n_discrep <- nrow(x$discrepancies)
  pct_discrep <- if (n_loci > 0) 100 * n_discrep / n_loci else 0
  cat(sprintf("Discrepant loci: %d (%.1f%% of panel loci)\n", n_discrep, pct_discrep))

  if (n_discrep > 0) {
    type_counts <- table(x$discrepancies$discrepancy_type)
    for (tn in names(type_counts)) {
      cat(sprintf("  - %s: %d\n", tn, type_counts[[tn]]))
    }

    show_n <- min(n_discrep, 5)
    cat(rule, "\n", sep = "")
    if (n_discrep == 1) {
      cat("Discrepant culprit locus:\n")
    } else {
      cat(sprintf("Top %d culprit loci with discrepancies:\n", show_n))
    }
    sub_disc <- x$discrepancies[seq_len(show_n), c("locus", "n_typed", "n_distinct_genotypes", "distinct_genotypes", "discrepancy_type")]
    print(as.data.frame(sub_disc), row.names = FALSE)
    if (n_discrep > show_n) {
      cat(sprintf("... and %d more discrepant loci in $discrepancies\n", n_discrep - show_n))
    }
  }

  cat(rule, "\n", sep = "")
  cat("Components: $summary (pairwise table), $discrepancies (culprit loci),\n")
  cat("            $distance_matrix (dist object), $ordination (PCoA coordinates).\n")
  cat("Use plot(x) to visualize replicate ordination.\n")
  invisible(x)
}

#' @export
print.replicate_trends <- function(x, ...) {
  g <- x$global
  cat(sprintf(
    "<replicate_trends> %d replicate set(s) (%d total samples) across %d loci\n",
    g$n_sets, g$n_samples_total, g$n_loci_total
  ))
  rule <- strrep("-", 50)
  cat(rule, "\n", sep = "")

  cat(sprintf(
    "Global concordance: mean %.1f%%, median %.1f%% (%d pairwise comparisons)\n",
    g$mean_concordance, g$median_concordance, g$n_pairs_total
  ))
  cat(sprintf(
    "Panel robustness:   %d of %d loci (%.1f%%) perfectly concordant across all sets\n",
    g$n_loci_clean, g$n_loci_total, g$pct_loci_clean
  ))

  # Outlier sets
  if (g$n_sets_flagged > 0) {
    cat(rule, "\n", sep = "")
    cat(sprintf(
      "! %d replicate set(s) flagged below %g%% concordance (potential mix-up):\n",
      g$n_sets_flagged, x$config$min_concordance
    ))
    bad_sets <- x$set_summary[x$set_summary$flagged_outlier, ]
    for (i in seq_len(nrow(bad_sets))) {
      cat(sprintf(
        "    * %s: min match %.1f%%, mean match %.1f%% (%d discrepant loci)\n",
        bad_sets$set_name[i], bad_sets$min_pct_match[i],
        bad_sets$mean_pct_match[i], bad_sets$n_discrepant_loci[i]
      ))
    }
  }

  # Chronic culprit loci
  flagged_loci <- x$locus_trends[x$locus_trends$flagged, ]
  cat(rule, "\n", sep = "")
  if (nrow(flagged_loci) > 0) {
    cat(sprintf(
      "Chronic culprit loci (discrepant in >= %.0f%% of sets): %d flagged\n",
      100 * x$config$locus_error_threshold, nrow(flagged_loci)
    ))
    show_n <- min(nrow(flagged_loci), 5)
    sub_l <- flagged_loci[seq_len(show_n), c("locus", "n_sets_discrepant", "n_sets_evaluated", "discrepancy_rate", "n_dropout", "n_discordant")]
    names(sub_l)[names(sub_l) == "discrepancy_rate"] <- "discrep_rate"
    sub_l$discrep_rate <- sprintf("%.1f%%", sub_l$discrep_rate)
    print(as.data.frame(sub_l), row.names = FALSE)
    if (nrow(flagged_loci) > show_n) {
      cat(sprintf("... and %d more flagged loci in $locus_trends\n", nrow(flagged_loci) - show_n))
    }
  } else {
    disc_any <- x$locus_trends[x$locus_trends$n_sets_discrepant > 0, ]
    if (nrow(disc_any) > 0) {
      cat(sprintf(
        "No chronic culprit loci flagged (none exceeded %.0f%% of sets). %d loci had occasional discrepancies.\n",
        100 * x$config$locus_error_threshold, nrow(disc_any)
      ))
    } else {
      cat("No discrepant loci found across any replicate sets.\n")
    }
  }

  cat(rule, "\n", sep = "")
  cat("Components: $locus_trends (all loci), $set_summary (per-set stats),\n")
  cat("            $global (panel metrics), $all_discrepancies (raw discrepancies).\n")
  cat("Use plot(x, type = 'culprits') or plot(x, type = 'concordance') to visualize.\n")
  invisible(x)
}


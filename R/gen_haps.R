# GEN_HAPS process port.
#
# Reimplements the chinook_gt Nextflow pipeline's GEN_HAPS process
# (modules/microhaplot.nf -> bin/gen_haps.R, which calls the microhaplotopia
# package) as dependency-free base R, so genoCompaR users can regenerate or
# re-parameterize a genotype sheet from raw haplotype calls without installing
# microhaplotopia or tidyverse. Every step below was checked against the real
# pipeline (microhaplotopia@emc-edits) on real data before being ported; see
# the CLAUDE.md "gen_haps()" section for the verified quirks this
# deliberately preserves rather than "fixes".

.gen_haps_required_cols <- c("group", "indiv.ID", "locus", "haplo", "depth",
                              "allele.balance", "rank")

# Per-locus/per-sample summary stats over a long-format call table, base R
# equivalent of microhaplotopia::summarize_data(). group_var is "locus" or
# "indiv.ID".
.gen_haps_summarize <- function(df, group_var) {
  if (nrow(df) == 0) {
    return(data.frame(
      key = character(0), n_samples = integer(0), n_loci = integer(0),
      n_haplotypes = integer(0), group_ids = character(0),
      mean_depth = numeric(0), median_depth = numeric(0),
      range_haplodepth = character(0), stringsAsFactors = FALSE
    ))
  }
  idx_by_group <- split(seq_len(nrow(df)), df[[group_var]])
  rows <- lapply(names(idx_by_group), function(g) {
    idx <- idx_by_group[[g]]
    d <- df$depth[idx]
    data.frame(
      key = g,
      n_samples = length(unique(df$indiv.ID[idx])),
      n_loci = length(unique(df$locus[idx])),
      n_haplotypes = length(unique(df$haplo[idx])),
      group_ids = paste(unique(df$group[idx]), collapse = ", "),
      mean_depth = mean(d),
      median_depth = stats::median(d),
      range_haplodepth = paste(range(d), collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Row count per locus (not distinct-sample count) among nxa - mirrors
# nxa %>% group_by(locus) %>% summarise(count = n()) %>% arrange(desc(count)).
.gen_haps_count_rows <- function(values) {
  if (length(values) == 0) {
    return(data.frame(key = character(0), count = integer(0), stringsAsFactors = FALSE))
  }
  tab <- table(values)
  out <- data.frame(key = names(tab), count = as.integer(tab), stringsAsFactors = FALSE)
  out[order(-out$count), ]
}

# Distinct-count of `value` per level of `by` - mirrors
# xtralleles %>% group_by(by) %>% summarise(count = n_distinct(value)).
.gen_haps_count_distinct <- function(value, by) {
  if (length(by) == 0) {
    return(data.frame(key = character(0), count = integer(0), stringsAsFactors = FALSE))
  }
  tab <- tapply(value, by, function(v) length(unique(v)))
  data.frame(key = names(tab), count = as.integer(tab), stringsAsFactors = FALSE)
}

#' Convert raw unfiltered haplotype calls into wide genotype layout
#'
#' Ports the `GEN_HAPS` process of the `chinook_gt` Nextflow pipeline
#' (`modules/microhaplot.nf` / `bin/gen_haps.R`, which calls the
#' `microhaplotopia` package) into a dependency-free R function: filters raw
#' per-allele haplotype calls by depth and allele balance, drops calls with an
#' N/X allele or more than two surviving alleles, resolves duplicate sample
#' re-runs, fills in the second allele for homozygous calls, and pivots the
#' result into the same wide `<locus>.1`/`<locus>.2` layout [read_genotypes()]
#' and [as_geno_data()] use.
#'
#' @param raw_haps Either a path to an `*_observed_unfiltered_haplotype.csv`
#'   file (long format: one row per allele call), or an already-loaded data
#'   frame in that format. Must have columns `group`, `indiv.ID`, `locus`,
#'   `haplo`, `depth`, `allele.balance`, `rank`; any other columns (e.g.
#'   `sum.Phred.C`, `max.Phred.C`, `source`) are ignored.
#' @param haplotype_depth Minimum read depth for a single allele call to be
#'   kept. Defaults to `4`, matching the Nextflow pipeline's
#'   `params.haplotype_depth`.
#' @param total_depth Minimum summed depth across a sample's surviving calls
#'   at a locus for that locus to be kept. Defaults to `8`, matching
#'   `params.total_depth`.
#' @param allele_balance Minimum allele balance for a call to be kept.
#'   Defaults to `0.35`, matching `params.allele_balance`.
#' @param preserve_group Keep each sample's own group instead of overwriting
#'   it with `"ND"`. `FALSE` by default, matching the Nextflow pipeline
#'   exactly (see Details). When `TRUE`, a re-run sample takes the group that
#'   won duplicate resolution (Step 6), and a sample re-added in Step 8 (no
#'   calls survived filtering) keeps the group its first raw call carried.
#' @param verbose Print the step-by-step report as the run proceeds. The
#'   report is retained either way; see [gen_haps_report()].
#'
#' @return A `geno_haps` object: a [tibble::tibble] with `group`, `indiv.ID`,
#'   then `<locus>.1`/`<locus>.2` columns for every locus with at least one
#'   surviving call, carrying the QC report as a `"gen_haps_report"`
#'   attribute. `group` is `"ND"` on every row (see Details). Sample IDs are
#'   unique, so `as_geno_data(out)` (its defaults already match) moves the
#'   result into the rest of genoCompaR.
#'
#' @details
#' The run proceeds in the same stages as `gen_haps.R`: filter by depth and
#' allele balance; drop N/X calls; drop loci with more than two surviving
#' alleles ("contaminated"); resolve duplicate sample re-runs by keeping
#' whichever run called the most loci; fill in the second allele for
#' homozygous calls; pivot to wide format; and add back, as all-`NA` rows,
#' any sample filtered out of the data entirely.
#'
#' Two behaviors are ported exactly because the already-shipped output of this
#' pipeline depends on them, not because they are obviously correct:
#' - **An N/X call wipes the whole `(sample, locus)` pair, not just the
#'   offending allele.** If one of a heterozygote's two calls is `N` or `X`,
#'   the surviving real call is also dropped, and the locus is missing (`NA`)
#'   for that sample rather than reduced to a homozygous-looking single call.
#' - **`group` is set to `"ND"` on every output row by default**, regardless
#'   of the sample's original group. Verified against the pipeline's shipped
#'   output (every row is `"ND"`) - not a placeholder this function
#'   introduces. Set `preserve_group = TRUE` to keep each sample's own group
#'   instead.
#'
#' Duplicate resolution is always the pipeline's "keep the run with the most
#' loci called, drop the rest" behavior; it is not exposed as an argument
#' because the Nextflow process does not parameterize it either.
#'
#' @seealso [gen_haps_report()] for the QC report, [as_geno_data()] to move
#'   the result into the rest of genoCompaR, [paralog_scan()] to screen the
#'   result for collapsed paralogs.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "synthetic", "synthetic_unfiltered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' out <- gen_haps(path, verbose = FALSE)
#' out
#' gen_haps_report(out)$anomalies
#'
#' geno <- as_geno_data(out)
#' get_loci(geno)
#'
#' @export
gen_haps <- function(raw_haps,
                      haplotype_depth = 4,
                      total_depth = 8,
                      allele_balance = 0.35,
                      preserve_group = FALSE,
                      verbose = TRUE) {

  if (allele_balance < 0 || allele_balance > 1) {
    stop("`allele_balance` must be between 0 and 1 (got ", allele_balance, ").",
         call. = FALSE)
  }
  if (haplotype_depth < 0) {
    stop("`haplotype_depth` must be non-negative (got ", haplotype_depth, ").",
         call. = FALSE)
  }
  if (total_depth < 0) {
    stop("`total_depth` must be non-negative (got ", total_depth, ").",
         call. = FALSE)
  }

  # See conform_genotypes() in R/conform_genotypes.R for the rationale behind
  # this closure-based report collector: per-run state, not package state, so
  # two runs never see each other's anomalies, and verbose = FALSE loses
  # nothing but the console printing.
  lines <- character(0)
  anomalies <- character(0)

  emit <- function(fmt, ...) {
    msg <- if (...length() == 0L) fmt else sprintf(fmt, ...)
    lines <<- c(lines, msg)
    if (verbose) cat(msg, "\n", sep = "")
    invisible(NULL)
  }
  note <- function(fmt, ...) {
    msg <- sprintf(fmt, ...)
    anomalies <<- c(anomalies, msg)
    emit("[ANOMALY] %s", msg)
  }

  # --- Step 1: load & validate ----------------------------------------------
  emit("--- Step 1: Loading raw haplotype calls ---")
  if (is.character(raw_haps) && length(raw_haps) == 1L) {
    if (!file.exists(raw_haps)) {
      stop("Required input not found: ", raw_haps, call. = FALSE)
    }
    raw <- .read_chr(raw_haps)
  } else if (is.data.frame(raw_haps)) {
    raw <- as.data.frame(raw_haps, stringsAsFactors = FALSE)
    for (col in intersect(c("group", "indiv.ID", "locus", "haplo"), names(raw))) {
      raw[[col]] <- as.character(raw[[col]])
    }
  } else {
    stop("`raw_haps` must be a file path or a data frame.", call. = FALSE)
  }

  missing_cols <- setdiff(.gen_haps_required_cols, names(raw))
  if (length(missing_cols) > 0) {
    stop("`raw_haps` is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (nrow(raw) == 0) stop("`raw_haps` has no rows.", call. = FALSE)

  raw$depth <- as.numeric(raw$depth)
  raw$allele.balance <- as.numeric(raw$allele.balance)
  raw$rank <- as.numeric(raw$rank)

  emit("[INFO] %d allele calls | %d samples | %d loci",
       nrow(raw), length(unique(raw$indiv.ID)), length(unique(raw$locus)))
  emit("[INFO] Thresholds: haplotype_depth >= %s, total_depth >= %s, allele_balance >= %s",
       haplotype_depth, total_depth, allele_balance)

  # --- Step 2: depth / allele-balance filter --------------------------------
  emit("--- Step 2: Depth and allele-balance filter ---")
  keep <- raw$allele.balance >= allele_balance & raw$depth >= haplotype_depth
  hap_fil <- raw[keep, , drop = FALSE]

  key <- paste(hap_fil$group, hap_fil$indiv.ID, hap_fil$locus, sep = "\r")
  hap_fil$TotalDepth <- stats::ave(hap_fil$depth, key, FUN = sum)
  hap_fil$NAlleles <- stats::ave(seq_len(nrow(hap_fil)), key, FUN = length)
  hap_fil <- hap_fil[hap_fil$TotalDepth >= total_depth, , drop = FALSE]

  emit("[INFO] %d of %d calls survive depth/balance/total-depth filtering.",
       nrow(hap_fil), nrow(raw))

  # --- Step 3: N/X allele calls ---------------------------------------------
  # Verified against the real pipeline: this removes every row sharing the
  # (group, indiv.ID, locus) combo, not just the offending N/X row - so a
  # heterozygote with one good call and one N call loses both.
  emit("--- Step 3: N/X allele calls ---")
  nxa <- hap_fil[grepl("N|X", hap_fil$haplo), , drop = FALSE]
  nxa_count <- .gen_haps_count_rows(nxa$locus)
  names(nxa_count) <- c("locus", "count")
  emit("[INFO] %d call(s) carry an N/X allele, affecting %d locus/loci.",
       nrow(nxa), length(unique(nxa$locus)))

  nxa_key <- unique(paste(nxa$group, nxa$indiv.ID, nxa$locus, sep = "\r"))
  hap_key <- paste(hap_fil$group, hap_fil$indiv.ID, hap_fil$locus, sep = "\r")
  hap_fil_nxa <- hap_fil[!(hap_key %in% nxa_key), , drop = FALSE]
  emit("[INFO] %d sample-locus pair(s) dropped entirely (N/X wipes both alleles).",
       length(nxa_key))

  # --- Step 4: contaminated (>2-allele) loci --------------------------------
  # NAlleles is exact here (not stale): Step 3 removes a combo's rows
  # all-or-nothing, so a surviving row's NAlleles still matches its row count.
  emit("--- Step 4: Contaminated (>2-allele) loci ---")
  xtralleles <- hap_fil_nxa[hap_fil_nxa$NAlleles > 2, , drop = FALSE]
  xtralleles <- xtralleles[order(xtralleles$group, xtralleles$indiv.ID, xtralleles$locus), ,
                           drop = FALSE]

  xtraalleles_individuals <- .gen_haps_count_distinct(xtralleles$locus, xtralleles$indiv.ID)
  names(xtraalleles_individuals) <- c("indiv.ID", "count")
  xtraalleles_locus <- .gen_haps_count_distinct(xtralleles$indiv.ID, xtralleles$locus)
  names(xtraalleles_locus) <- c("locus", "count")
  xtraalleles_locus <- xtraalleles_locus[order(-xtraalleles_locus$count), , drop = FALSE]

  hap_fil1 <- hap_fil_nxa[hap_fil_nxa$NAlleles <= 2, , drop = FALSE]
  n_contaminated <- length(unique(paste(xtralleles$indiv.ID, xtralleles$locus, sep = "\r")))
  emit("[INFO] %d sample-locus pair(s) dropped for carrying more than 2 alleles.",
       n_contaminated)
  if (n_contaminated > 0) {
    note("%d sample-locus pair(s) carried more than 2 alleles and were dropped.",
         n_contaminated)
  }

  # --- Step 5: per-locus / per-sample depth diagnostics ---------------------
  emit("--- Step 5: Locus and sample depth summaries ---")
  locus_depth <- .gen_haps_summarize(hap_fil1, "locus")
  names(locus_depth)[1] <- "locus"
  locus_depth <- locus_depth[order(locus_depth$n_samples), , drop = FALSE]

  individual_depth <- .gen_haps_summarize(hap_fil1, "indiv.ID")
  names(individual_depth)[1] <- "indiv.ID"
  individual_depth <- individual_depth[order(individual_depth$mean_depth), , drop = FALSE]

  # --- Step 6: duplicate sample re-runs (resolve = "drop") ------------------
  emit("--- Step 6: Duplicate sample re-runs ---")
  distinct_pairs <- unique(hap_fil1[, c("group", "indiv.ID")])
  id_counts <- table(distinct_pairs$indiv.ID)
  dup_ids <- names(id_counts)[id_counts > 1]

  duplicate_samples <- distinct_pairs[distinct_pairs$indiv.ID %in% dup_ids, c("indiv.ID", "group"),
                                       drop = FALSE]
  duplicate_samples <- duplicate_samples[order(duplicate_samples$indiv.ID), , drop = FALSE]

  if (length(dup_ids) > 0) {
    note("%d sample ID(s) were re-run under more than one group tag; kept the run with the most loci called.",
         length(dup_ids))
    keep <- rep(TRUE, nrow(hap_fil1))
    for (id in dup_ids) {
      rws <- which(hap_fil1$indiv.ID == id)
      # Tie-break: alphabetically-first group wins, matching the real
      # pipeline's arrange(desc(N_Loci)) over group-then-indiv.ID-sorted rows.
      grps <- sort(unique(hap_fil1$group[rws]))
      n_loci_by_group <- vapply(grps, function(g) {
        length(unique(hap_fil1$locus[rws[hap_fil1$group[rws] == g]]))
      }, integer(1))
      winner <- grps[which.max(n_loci_by_group)]
      drop_idx <- rws[hap_fil1$group[rws] != winner]
      keep[drop_idx] <- FALSE
      emit("  %-14s kept group '%s' (%d loci); dropped group(s): %s",
           id, winner, max(n_loci_by_group),
           paste(setdiff(grps, winner), collapse = ", "))
    }
    hap_fil1 <- hap_fil1[keep, , drop = FALSE]
  } else {
    emit("[OK] No duplicate sample IDs present.")
  }

  # --- Step 7: fill homozygotes and pivot to wide layout --------------------
  emit("--- Step 7: Filling homozygous calls and pivoting to wide layout ---")
  ord <- order(hap_fil1$indiv.ID, hap_fil1$locus, hap_fil1$rank)
  hf <- hap_fil1[ord, , drop = FALSE]
  cell_key <- paste(hf$indiv.ID, hf$locus, sep = "\r")
  first_idx <- !duplicated(cell_key)
  # NAlleles <= 2 guarantees at most 2 rows per (indiv.ID, locus), so the
  # second occurrence (if any) is exactly the other allele.
  second_idx <- duplicated(cell_key)

  ids_sorted <- sort(unique(hf$indiv.ID))
  loci_sorted <- sort(unique(hf$locus))

  a1 <- hf[first_idx, , drop = FALSE]
  a2 <- hf[second_idx, , drop = FALSE]

  mat1 <- matrix(NA_character_, nrow = length(ids_sorted), ncol = length(loci_sorted),
                 dimnames = list(ids_sorted, loci_sorted))
  mat2 <- mat1
  if (nrow(a1) > 0) {
    mat1[cbind(match(a1$indiv.ID, ids_sorted), match(a1$locus, loci_sorted))] <- a1$haplo
  }
  if (nrow(a2) > 0) {
    mat2[cbind(match(a2$indiv.ID, ids_sorted), match(a2$locus, loci_sorted))] <- a2$haplo
  }
  hom <- is.na(mat2) & !is.na(mat1)
  mat2[hom] <- mat1[hom]

  locus_cols <- as.vector(rbind(paste0(loci_sorted, ".1"), paste0(loci_sorted, ".2")))
  wide <- matrix(NA_character_, nrow = length(ids_sorted), ncol = length(locus_cols))
  if (length(loci_sorted) > 0) {
    wide[, seq(1, ncol(wide), by = 2)] <- mat1
    wide[, seq(2, ncol(wide), by = 2)] <- mat2
  }
  wide_df <- as.data.frame(wide, stringsAsFactors = FALSE)
  names(wide_df) <- locus_cols
  wide_df <- cbind(data.frame(indiv.ID = ids_sorted, stringsAsFactors = FALSE), wide_df)

  emit("[INFO] %d samples x %d loci retained after filtering.",
       length(ids_sorted), length(loci_sorted))

  # --- Step 8: re-add samples filtered out entirely --------------------------
  emit("--- Step 8: Re-adding filtered-out samples ---")
  raw_ids <- unique(raw$indiv.ID)
  missing_ids <- raw_ids[!raw_ids %in% ids_sorted]
  missing_samples <- data.frame(
    indiv.ID = missing_ids,
    group = raw$group[match(missing_ids, raw$indiv.ID)],
    stringsAsFactors = FALSE
  )
  emit("[INFO] %d sample(s) had no calls survive filtering; re-added as all-NA rows.",
       length(missing_ids))

  if (length(missing_ids) > 0) {
    missing_block <- as.data.frame(
      matrix(NA_character_, nrow = length(missing_ids), ncol = length(locus_cols)),
      stringsAsFactors = FALSE
    )
    names(missing_block) <- locus_cols
    missing_block <- cbind(
      data.frame(indiv.ID = missing_ids, stringsAsFactors = FALSE), missing_block
    )
    wide_df <- rbind(wide_df, missing_block)
  }

  # --- Step 9: finalize -------------------------------------------------------
  emit("--- Step 9: Finalizing ---")
  if (preserve_group) {
    # Each surviving indiv.ID has exactly one group left in hf after Step 6's
    # duplicate resolution, so any of its rows is representative.
    survived_group <- hf$group[match(ids_sorted, hf$indiv.ID)]
    wide_df$group <- c(survived_group, missing_samples$group)
    emit("[INFO] preserve_group = TRUE: output carries each sample's own group.")
  } else {
    # Matches the pipeline exactly: group is "ND" on every row, not just the
    # re-added ones. See @details above.
    wide_df$group <- "ND"
    emit("[INFO] preserve_group = FALSE (default): group is 'ND' on every row, matching the pipeline.")
  }
  wide_df <- wide_df[, c("group", "indiv.ID", locus_cols), drop = FALSE]
  rownames(wide_df) <- NULL

  emit("[INFO] Final table: %d samples x %d loci.", nrow(wide_df), length(loci_sorted))
  dropped_loci <- setdiff(unique(raw$locus), loci_sorted)
  if (length(dropped_loci) > 0) {
    note("%d locus/loci had no calls survive filtering and are absent from the output: %s",
         length(dropped_loci),
         paste(utils::head(sort(dropped_loci), 10), collapse = ", "))
  }
  if (nrow(wide_df) > 0 && length(missing_ids) / nrow(wide_df) > 0.2) {
    note("%.1f%% of samples had no calls survive filtering.",
         100 * length(missing_ids) / nrow(wide_df))
  }

  emit("")
  if (length(anomalies) == 0) {
    emit("[OK] No anomalies noted.")
  } else {
    emit("%d condition(s) noted above may need a human look.", length(anomalies))
  }

  report <- structure(
    list(
      text = lines,
      anomalies = anomalies,
      params = list(haplotype_depth = haplotype_depth, total_depth = total_depth,
                    allele_balance = allele_balance, preserve_group = preserve_group),
      nxa = tibble::as_tibble(nxa),
      nxa_count = tibble::as_tibble(nxa_count),
      xtraalleles = tibble::as_tibble(xtralleles),
      xtraalleles_individuals = tibble::as_tibble(xtraalleles_individuals),
      xtraalleles_locus = tibble::as_tibble(xtraalleles_locus),
      locus_depth = tibble::as_tibble(locus_depth),
      individual_depth = tibble::as_tibble(individual_depth),
      duplicate_samples = tibble::as_tibble(duplicate_samples),
      missing_samples = tibble::as_tibble(missing_samples)
    ),
    class = "gen_haps_report"
  )

  out <- tibble::as_tibble(wide_df)
  structure(
    out,
    class = c("geno_haps", class(out)),
    gen_haps_report = report
  )
}

#' The QC report from a gen_haps() run
#'
#' @param x A `geno_haps` object from [gen_haps()].
#' @return A `gen_haps_report`: a list of `text` (every line of the report, in
#'   order), `anomalies` (conditions worth a human look), `params` (the
#'   thresholds used), and the same diagnostic tables `gen_haps.R` writes to
#'   disk: `nxa`, `nxa_count`, `xtraalleles`, `xtraalleles_individuals`,
#'   `xtraalleles_locus`, `locus_depth`, `individual_depth`,
#'   `duplicate_samples`, and `missing_samples`. Printing it re-renders the
#'   report.
#' @seealso [gen_haps()]
#' @examples
#' path <- system.file(
#'   "extdata", "synthetic", "synthetic_unfiltered_haplotypes.csv",
#'   package = "genoCompaR"
#' )
#' out <- gen_haps(path, verbose = FALSE)
#' gen_haps_report(out)$missing_samples
#'
#' @export
gen_haps_report <- function(x) {
  report <- attr(x, "gen_haps_report")
  if (is.null(report)) {
    stop("`x` must be a geno_haps object (see gen_haps()).", call. = FALSE)
  }
  report
}

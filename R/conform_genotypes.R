# CFS -> DWR haplotype conversion.
#
# The QC here is driven by the reference VCFs rather than by lengths inferred
# from the reference CSV: see R/vcf_alleles.R for why that distinction is the
# whole point of the pipeline.

#' Translate a CFS marker column name into its DWR equivalent.
#'
#' @param name Character vector of CFS marker column names.
#' @return Character vector of DWR column names.
#' @noRd
.transform_name <- function(name) {
  vapply(name, function(nm) {
    # a. Suffix shift: CFS uses "locus" / "locus.1" for the two alleles,
    #    DWR uses "locus.1" / "locus.2".
    if (grepl("\\.1$", nm)) {
      base <- sub("\\.1$", "", nm)
      suffix <- ".2"
    } else {
      base <- nm
      suffix <- ".1"
    }

    # b. NC_ genomic coordinates
    #    NC_037099.1.62937268.62937373 -> NC_037099.1:62937268-62937373
    if (grepl("^NC_", base)) {
      base <- sub("^NC_(\\d+\\.\\d+)\\.(\\d+)\\.(\\d+)$", "NC_\\1:\\2-\\3", base)
    }
    # c. Ots_ / OkiOts_ dash formatting
    else if (grepl("^Ots_|^OkiOts_", base)) {
      if (grepl("u07", base)) {
        # u07 keeps its trailing dot: Ots_u07.07.161 -> Ots_u07-07.161
        base <- sub("\\.", "-", base)
      } else {
        # Everything else dashes throughout:
        #   Ots_111312.435    -> Ots_111312-435
        #   Ots_BMP.2.SNP1    -> Ots_BMP-2-SNP1
        #   OkiOts_120255.113 -> OkiOts_120255-113
        base <- gsub("\\.", "-", base)
      }
    }

    paste0(base, suffix)
  }, character(1), USE.NAMES = FALSE)
}

#' Convert a CFS genotype export into the DWR haplotype layout
#'
#' Reads a CFS-format Chinook genotype export, harmonizes its marker names to
#' the DWR convention, validates every allele against reference VCFs, resolves
#' repeated sample IDs, and returns a table conforming exactly to the column
#' layout of a DWR reference file - the format downstream tools such as
#' CKMRsim and CKMR-shiny expect.
#'
#' @param cfs_file Path to the CFS genotype export (CSV). One row per sample;
#'   marker alleles named `<locus>` and `<locus>.1`.
#' @param dwr_file Path to a DWR reference genotype file (CSV). Supplies the
#'   exact output column layout and the observed-allele catalogue. A file
#'   leading with `group`/`indiv.ID` instead of `SAMPLE_ID` is accepted;
#'   `group` is dropped and `indiv.ID` renamed.
#' @param vcf_files Character vector of reference VCF paths. A path that does
#'   not exist is reported and skipped rather than raising, but shrinks QC
#'   coverage.
#' @param id_col Name of the CFS column holding sample identifiers. Defaults
#'   to `"indiv"`.
#' @param meta_cols Non-genotype columns in the CFS export. Only `id_col` is
#'   required; the rest are dropped if present and noted if absent, so a
#'   changed metadata block does not stop the run.
#' @param repair_offsets Move genotype blocks the offset scan diagnosed as
#'   shifted back to the loci they belong to. `FALSE` by default: repair
#'   rewrites data, so the intended workflow is to run once with it off, read
#'   the offset report, and turn it on only once the diagnosis is accepted.
#' @param null_invalid_alleles Blank alleles the VCFs say are impossible at
#'   their locus and that no shift explains. `FALSE` by default - a real allele
#'   absent from the reference panel looks the same as an artefact.
#' @param resolve_duplicate_ids Collapse repeated sample IDs (CFS re-runs) to
#'   one row, keeping whichever row carries the most called genotypes. IDs
#'   whose every row is empty are dropped. `TRUE` by default.
#' @param locus_shifts Offsets, in whole loci, the scan searches for. A
#'   locus's two allele columns move together, so shifts are never half a
#'   locus.
#' @param min_offset_run Number of consecutive loci that must share one offset
#'   before the run is reported.
#' @param low_call_rate Call rate below which a sample is called out in the
#'   report.
#' @param max_het_excess,min_het_samples Passed to [paralog_scan()] as
#'   `threshold` and `min_n` for the collapsed-paralog check.
#' @param max_replicate_discordance Fraction of shared calls at which two runs
#'   of one sample are considered to disagree, raising an anomaly.
#' @param verbose Print the QC report as the run proceeds. The report is
#'   retained either way; see [conform_report()].
#'
#' @return A `geno_conformed` object: a [tibble::tibble] with the DWR
#'   reference's exact columns, in order, carrying the QC report as a
#'   `"conform_report"` attribute.
#'
#' @details
#' The run proceeds in eight steps, each announced in the report: load; check
#' allele-column pairing; harmonize marker names; scan for column offsets;
#' optionally repair them; null alleles that fail validity or width; conform to
#' the DWR template and resolve duplicate IDs; and report.
#'
#' Allele validity comes from the VCFs rather than from column widths because a
#' CFS export's columns can be shifted one locus off from their header with no
#' visible symptom, if the shifted-in locus happens to have the same haplotype
#' length as its neighbour. Checking each base against the VCF's permitted
#' bases at that position catches the shift anyway.
#'
#' Anything the run notices but has no specific rule for is collected as an
#' anomaly and listed at the end, so a novel fault surfaces instead of passing
#' silently.
#'
#' The result carries duplicate sample IDs when `resolve_duplicate_ids` is
#' `FALSE`, so it is a plain table rather than a `geno_data` object. Convert it
#' with `as_geno_data(out, id_col = "SAMPLE_ID", group_col = NULL)` to compare
#' samples with the rest of the package.
#'
#' @seealso [conform_report()] for the QC report, [paralog_scan()] for the
#'   heterozygote-excess check run as part of it, and [as_geno_data()] to move
#'   the result into the rest of genoCompaR.
#'
#' @examples
#' syn <- function(f) {
#'   system.file("extdata", "synthetic", f, package = "genoCompaR")
#' }
#' out <- conform_genotypes(
#'   syn("synthetic_cfs_clean.csv"),
#'   syn("synthetic_dwr_reference.csv"),
#'   syn(c("synthetic_full.vcf", "synthetic_panel.vcf")),
#'   verbose = FALSE
#' )
#' out
#' conform_report(out)$anomalies
#'
#' geno <- as_geno_data(out, id_col = "SAMPLE_ID", group_col = NULL)
#' length(get_loci(geno))
#'
#' @export
conform_genotypes <- function(cfs_file,
                               dwr_file,
                               vcf_files,
                               id_col = "indiv",
                               meta_cols = c("indiv", "repunit", "collection",
                                             "sample_type", "rosa_genotypes_str"),
                               repair_offsets = FALSE,
                               null_invalid_alleles = FALSE,
                               resolve_duplicate_ids = TRUE,
                               locus_shifts = c(-1L, 1L, -2L, 2L),
                               min_offset_run = 3L,
                               low_call_rate = 0.5,
                               max_het_excess = 0.25,
                               min_het_samples = 50L,
                               max_replicate_discordance = 0.01,
                               verbose = TRUE) {

  # The report is accumulated as text and returned, so a caller who ran with
  # verbose = FALSE (or in a script, or in a test) can still read every line
  # the run produced. Anomalies are collected per run in this closure rather
  # than in a package-level variable, so two runs cannot see each other's.
  lines <- character(0)
  anomalies <- character(0)

  emit <- function(fmt, ...) {
    msg <- if (...length() == 0L) fmt else sprintf(fmt, ...)
    lines <<- c(lines, msg)
    if (verbose) cat(msg, "\n", sep = "")
    invisible(NULL)
  }
  emit_obj <- function(obj, ...) {
    txt <- utils::capture.output(print(obj, ...))
    lines <<- c(lines, txt)
    # One trailing newline, not two: the console and `lines` must agree, so
    # that a run with verbose = FALSE loses nothing but the printing.
    if (verbose) cat(paste0(paste(txt, collapse = "\n"), "\n"))
    invisible(NULL)
  }
  note <- function(fmt, ...) {
    msg <- sprintf(fmt, ...)
    anomalies <<- c(anomalies, msg)
    emit("[ANOMALY] %s", msg)
  }

  # --- Step 1: load ----------------------------------------------------------
  emit("--- Step 1: Loading files ---")
  for (p in c(cfs_file, dwr_file)) {
    if (!file.exists(p)) {
      stop("Required input not found: ", p,
           "\n  (working directory: ", getwd(), ")", call. = FALSE)
    }
  }
  # colClasses = "character" keeps all-numeric marker columns from becoming
  # integers, which would make nchar() report digit counts instead of
  # haplotype lengths and corrupt the width QC below.
  cfs_raw <- .read_chr(cfs_file)
  dwr_ref <- .read_chr(dwr_file)
  vcf_alleles <- .load_vcf_alleles(vcf_files, note = note, emit = emit)

  emit("[INFO] CFS: %d rows x %d cols | DWR: %d rows x %d cols | VCF: %d loci",
       nrow(cfs_raw), ncol(cfs_raw), nrow(dwr_ref), ncol(dwr_ref),
       length(vcf_alleles))

  if (nrow(cfs_raw) == 0) stop("CFS file has no data rows: ", cfs_file, call. = FALSE)

  # Some DWR exports lead with "group,indiv.ID,..." instead of "SAMPLE_ID,...".
  # group carries no genotype information, so it is dropped; indiv.ID is the
  # sample identifier under a different name, so it is renamed in place.
  if (!"SAMPLE_ID" %in% colnames(dwr_ref) &&
      length(colnames(dwr_ref)) >= 2 &&
      identical(colnames(dwr_ref)[1:2], c("group", "indiv.ID"))) {
    note(paste("DWR reference used 'group'/'indiv.ID' header instead of",
               "SAMPLE_ID -- group dropped, indiv.ID renamed to SAMPLE_ID."))
    dwr_ref <- dwr_ref[, setdiff(colnames(dwr_ref), "group"), drop = FALSE]
    colnames(dwr_ref)[colnames(dwr_ref) == "indiv.ID"] <- "SAMPLE_ID"
  }
  if (!"SAMPLE_ID" %in% colnames(dwr_ref)) {
    stop("DWR reference has no SAMPLE_ID column; cannot build the output template.",
         call. = FALSE)
  }
  dwr_cols <- colnames(dwr_ref)
  marker_cols_dwr <- setdiff(dwr_cols, "SAMPLE_ID")

  if (!id_col %in% colnames(cfs_raw)) {
    stop("CFS file has no '", id_col, "' column to use as SAMPLE_ID.",
         "\n  Columns seen: ",
         paste(utils::head(colnames(cfs_raw), 8), collapse = ", "), " ...",
         call. = FALSE)
  }
  # Absent metadata is survivable; only note it so a schema change is visible.
  missing_meta <- setdiff(meta_cols, colnames(cfs_raw))
  if (length(missing_meta) > 0) {
    note("CFS metadata columns absent (continuing): %s",
         paste(missing_meta, collapse = ", "))
  }

  cfs_markers_raw <- setdiff(colnames(cfs_raw), meta_cols)
  if (length(cfs_markers_raw) == 0) stop("CFS file has no marker columns.", call. = FALSE)

  blank_cols <- colnames(cfs_raw)[colnames(cfs_raw) == "" | is.na(colnames(cfs_raw))]
  if (length(blank_cols) > 0) {
    note("%d CFS columns have blank names -- possible trailing delimiters.",
         length(blank_cols))
  }
  n_samples <- nrow(cfs_raw)

  # --- Step 2: structural check on allele pairing ----------------------------
  # A malformed header shows up here before it can corrupt the suffix shift.
  emit("--- Step 2: Checking allele column pairing ---")
  seconds <- cfs_markers_raw[grepl("\\.1$", cfs_markers_raw)]
  firsts <- cfs_markers_raw[!grepl("\\.1$", cfs_markers_raw)]
  orphan_second <- seconds[!(sub("\\.1$", "", seconds) %in% firsts)]
  orphan_first <- firsts[!(paste0(firsts, ".1") %in% seconds)]

  if (length(cfs_markers_raw) %% 2 != 0) {
    note("Odd number of marker columns (%d) -- alleles do not pair evenly.",
         length(cfs_markers_raw))
  }
  if (length(orphan_second) > 0 || length(orphan_first) > 0) {
    # Unpaired allele columns mean the locus-level offset scan cannot trust its
    # column grouping, so this is a hard anomaly rather than a passing note.
    note("%d allele columns have no partner; locus grouping is unreliable.",
         length(orphan_second) + length(orphan_first))
    if (length(orphan_first) > 0) {
      emit("  First alleles with no .1 partner: %s",
           paste(utils::head(orphan_first, 5), collapse = ", "))
    }
    if (length(orphan_second) > 0) {
      emit("  Second alleles with no base:      %s",
           paste(utils::head(orphan_second, 5), collapse = ", "))
    }
  } else {
    emit("[OK] %d marker columns form %d complete allele pairs.",
         length(cfs_markers_raw), length(cfs_markers_raw) %/% 2)
  }

  # --- Step 3: rename --------------------------------------------------------
  emit("--- Step 3: Harmonizing CFS marker names ---")
  cfs_renamed <- cfs_raw[, cfs_markers_raw, drop = FALSE]
  colnames(cfs_renamed) <- .transform_name(cfs_markers_raw)

  dup_names <- unique(colnames(cfs_renamed)[duplicated(colnames(cfs_renamed))])
  if (length(dup_names) > 0) {
    emit("[WARNING] Renaming produced %d duplicate column names: %s",
         length(dup_names), paste(utils::head(dup_names, 5), collapse = ", "))
  }

  matches <- intersect(colnames(cfs_renamed), marker_cols_dwr)
  unmatched <- setdiff(colnames(cfs_renamed), marker_cols_dwr)

  emit("[INFO] %d of %d CFS markers mapped onto DWR names.",
       length(matches), ncol(cfs_renamed))
  if (length(unmatched) > 0) {
    # Do not let these pass quietly. A rename bug looks exactly like a marker
    # that legitimately is not in the DWR panel.
    emit("[WARNING] %d CFS markers did NOT map and will be dropped:",
         length(unmatched))
    emit_obj(utils::head(unmatched, 20))
    if (length(unmatched) > 20) {
      emit("  ... and %d more.", length(unmatched) - 20)
    }
  }

  # --- Step 4: locus-level offset detection ---------------------------------
  # Runs before any cleaning, on the untouched data. A locus's two allele
  # columns move together, so offsets are searched in whole loci: shift d means
  # the genotype sitting at locus i actually belongs to locus i+d.
  emit("--- Step 4: Column-offset scan (VCF-validated) ---")

  permitted <- .build_permitted(vcf_alleles, dwr_ref, emit = emit)

  loci <- sub("\\.1$", "", colnames(cfs_renamed)[c(TRUE, FALSE)])
  n_loci <- length(loci)
  has_vcf <- loci %in% names(permitted)
  emit("[INFO] %d of %d loci are described by the VCF.", sum(has_vcf), n_loci)

  # ok[i, L] = sample i's genotype at locus L is valid for locus L+d.
  validity_at <- function(d) {
    m <- matrix(NA, nrow = n_samples, ncol = n_loci)
    for (L in seq_len(n_loci)) {
      tgt <- L + d
      if (tgt < 1 || tgt > n_loci || !has_vcf[tgt]) next
      perm <- permitted[[loci[tgt]]]
      a1 <- cfs_renamed[[paste0(loci[L], ".1")]]
      a2 <- cfs_renamed[[paste0(loci[L], ".2")]]
      called <- .is_called(a1) | .is_called(a2)
      good <- (!.is_called(a1) | .valid_allele(a1, perm)) &
        (!.is_called(a2) | .valid_allele(a2, perm))
      m[, L] <- ifelse(called, good, NA)
    }
    m
  }

  ok_self <- validity_at(0L)
  ok_shift <- lapply(locus_shifts, validity_at)
  names(ok_shift) <- as.character(locus_shifts)

  # Per sample and locus, the offset that explains the genotype (0 = in place).
  shift_mat <- matrix(NA_integer_, nrow = n_samples, ncol = n_loci)
  shift_mat[which(ok_self %in% TRUE)] <- 0L
  for (d in locus_shifts) {
    take <- which(!(ok_self %in% TRUE) & (ok_shift[[as.character(d)]] %in% TRUE) &
                    is.na(shift_mat))
    shift_mat[take] <- d
  }

  # Collapse into runs of consecutive loci sharing one non-zero offset.
  offset_runs <- list()
  for (i in seq_len(n_samples)) {
    s <- shift_mat[i, ]
    s[is.na(s)] <- 0L   # uncalled loci do not break a run
    r <- rle(s)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    for (k in seq_along(r$values)) {
      if (r$values[k] == 0L || r$lengths[k] < min_offset_run) next
      offset_runs[[length(offset_runs) + 1]] <- data.frame(
        sample_row = i,
        SAMPLE_ID = cfs_raw[[id_col]][i],
        from_locus = loci[starts[k]],
        to_locus = loci[ends[k]],
        n_loci = r$lengths[k],
        shift = r$values[k],
        start_idx = starts[k],
        end_idx = ends[k],
        stringsAsFactors = FALSE
      )
    }
  }

  runs <- NULL
  offset_table <- NULL
  offset_samples <- character(0)
  if (length(offset_runs) > 0) {
    runs <- do.call(rbind, offset_runs)
    offset_samples <- unique(runs$SAMPLE_ID)

    emit("[WARNING] %d samples carry a shifted block (%.1f%% of samples).",
         length(offset_samples), 100 * length(offset_samples) / n_samples)

    # Run boundaries wander with each sample's missing calls, so the same
    # underlying corruption appears as many near-identical blocks. Summarise
    # per locus instead: that shows the affected region and its extent plainly.
    affected <- matrix(FALSE, nrow = n_samples, ncol = n_loci)
    for (k in seq_len(nrow(runs))) {
      affected[runs$sample_row[k], runs$start_idx[k]:runs$end_idx[k]] <- TRUE
    }
    n_aff <- colSums(affected)
    if (any(n_aff > 0)) {
      dom <- vapply(seq_len(n_loci), function(L) {
        if (n_aff[L] == 0) return(NA_integer_)
        s <- shift_mat[affected[, L], L]
        s <- s[!is.na(s) & s != 0L]
        if (length(s) == 0) return(NA_integer_)
        as.integer(names(sort(table(s), decreasing = TRUE))[1])
      }, integer(1))
      tab <- data.frame(locus = loci, samples = n_aff, shift = dom,
                        stringsAsFactors = FALSE)
      tab <- tab[tab$samples > 0, ]
      offset_table <- tab[order(-tab$samples), ]
      emit("Loci most often shifted (locus, samples affected, dominant shift):")
      emit_obj(utils::head(offset_table, 25), row.names = FALSE)
      if (nrow(offset_table) > 25) {
        emit("  ... and %d more loci.", nrow(offset_table) - 25)
      }
      emit("A contiguous run at one shift affecting most samples is a systematic")
      emit("offset in the CFS export; scattered loci at low counts are per-sample.")
    }
  } else {
    emit("[OK] No shifted blocks detected.")
  }

  # Genotypes the VCF rejects at their own locus and that no shift explains.
  n_unexplained <- sum(is.na(shift_mat) & !is.na(ok_self))
  emit("[INFO] %d sample-loci are invalid in place and unexplained by any shift.",
       n_unexplained)

  # --- Step 5: optional repair ----------------------------------------------
  if (repair_offsets && length(offset_runs) > 0) {
    emit("--- Step 5: Repairing shifted blocks ---")
    moved <- 0L
    for (k in seq_len(nrow(runs))) {
      i <- runs$sample_row[k]
      d <- runs$shift[k]
      src <- runs$start_idx[k]:runs$end_idx[k]
      dst <- src + d
      keep <- dst >= 1 & dst <= n_loci
      # Read the whole block first: source and destination ranges overlap.
      v1 <- vapply(loci[src[keep]], function(L) cfs_renamed[i, paste0(L, ".1")],
                   character(1))
      v2 <- vapply(loci[src[keep]], function(L) cfs_renamed[i, paste0(L, ".2")],
                   character(1))
      for (L in loci[src]) {
        cfs_renamed[i, paste0(L, ".1")] <- NA_character_
        cfs_renamed[i, paste0(L, ".2")] <- NA_character_
      }
      for (m in seq_along(dst[keep])) {
        cfs_renamed[i, paste0(loci[dst[keep][m]], ".1")] <- v1[m]
        cfs_renamed[i, paste0(loci[dst[keep][m]], ".2")] <- v2[m]
      }
      moved <- moved + sum(keep)
    }
    emit("[INFO] Moved %d sample-loci back into place across %d blocks.",
         moved, nrow(runs))
    emit("[INFO] Genotypes shifted off either end of the panel are unrecoverable.")
  } else if (length(offset_runs) > 0) {
    emit("--- Step 5: Repair skipped ---")
    emit("[INFO] Shifted blocks left in place; set repair_offsets = TRUE to move them.")
    emit("[INFO] Left unrepaired, these genotypes are wrong and mostly get nulled in Step 6.")
    note(paste("%d samples carry shifted genotype blocks and repair is OFF;",
               "their affected loci are wrong in the output."),
         length(offset_samples))
  }

  # --- Step 6: allele validity and length QC --------------------------------
  # CKMRsim requires every allele at a locus to be one fixed width. Anything the
  # VCF rejects outright, or that disagrees with the reference width, is nulled.
  emit("--- Step 6: Allele validity / length QC ---")

  nulled_invalid <- 0L
  if (null_invalid_alleles) {
    for (L in loci[has_vcf]) {
      perm <- permitted[[L]]
      for (col in paste0(L, c(".1", ".2"))) {
        if (!col %in% colnames(cfs_renamed)) next
        bad <- .is_called(cfs_renamed[[col]]) &
          !.valid_allele(cfs_renamed[[col]], perm)
        if (any(bad)) {
          cfs_renamed[bad, col] <- NA_character_
          nulled_invalid <- nulled_invalid + sum(bad)
        }
      }
    }
    emit("[INFO] Nulled %d alleles the VCF rejects (null_invalid_alleles = TRUE).",
         nulled_invalid)
  }

  # Expected width: VCF site count where known, else the reference CSV.
  true_length <- list()
  for (col in matches) {
    L <- sub("\\.[12]$", "", col)
    if (!is.null(vcf_alleles[[L]])) {
      true_length[[col]] <- length(vcf_alleles[[L]])
    } else {
      v <- dwr_ref[[col]]
      lens <- unique(nchar(v[.is_called(v)]))
      if (length(lens) == 1) true_length[[col]] <- lens
    }
  }

  nulled_per_col <- stats::setNames(integer(length(matches)), matches)
  called_per_col <- stats::setNames(integer(length(matches)), matches)

  for (col in matches) {
    vals <- cfs_renamed[[col]]
    called <- .is_called(vals)
    called_per_col[col] <- sum(called)
    tl <- true_length[[col]]
    if (is.null(tl)) next
    bad <- called & nchar(vals) != tl
    if (any(bad)) {
      cfs_renamed[[col]][bad] <- NA_character_
      nulled_per_col[col] <- sum(bad)
    }
  }

  emit("[INFO] Nulled %d alleles across %d columns on width mismatch.",
       sum(nulled_per_col), sum(nulled_per_col > 0))

  frac_nulled <- ifelse(called_per_col > 0, nulled_per_col / called_per_col, 0)
  gutted_loci <- sort(unique(sub("\\.[12]$", "", names(frac_nulled)[frac_nulled > 0.5])))

  # --- Step 7: conform to the DWR template ----------------------------------
  emit("--- Step 7: Conforming to DWR template ---")
  conformed_df <- as.data.frame(
    matrix(NA_character_, nrow = n_samples, ncol = length(dwr_cols)),
    stringsAsFactors = FALSE
  )
  colnames(conformed_df) <- dwr_cols
  conformed_df$SAMPLE_ID <- cfs_raw[[id_col]]
  conformed_df[, matches] <- cfs_renamed[, matches, drop = FALSE]

  blank_ids <- !.is_called(conformed_df$SAMPLE_ID)
  if (any(blank_ids)) {
    note("%d rows have a blank SAMPLE_ID.", sum(blank_ids))
  }

  # --- Step 7b: resolve duplicate SAMPLE_IDs --------------------------------
  # CFS re-runs a sample under its original ID, so an ID can appear more than
  # once: typically one failed (empty) row and one good one, but sometimes two
  # real genotypes that disagree. Keep the row carrying the most calls; drop
  # IDs whose every row is empty, since they contribute nothing either way.
  dropped_empty_ids <- character(0)
  if (resolve_duplicate_ids && any(duplicated(conformed_df$SAMPLE_ID))) {
    emit("--- Step 7b: Resolving duplicate SAMPLE_IDs ---")

    called_counts <- rowSums(!is.na(conformed_df[, marker_cols_dwr, drop = FALSE]))
    dup_ids <- unique(conformed_df$SAMPLE_ID[duplicated(conformed_df$SAMPLE_ID)])
    emit("[INFO] %d SAMPLE_IDs appear more than once (%d rows involved).",
         length(dup_ids), sum(conformed_df$SAMPLE_ID %in% dup_ids))

    keep <- rep(TRUE, nrow(conformed_df))
    for (id in dup_ids) {
      rws <- which(conformed_df$SAMPLE_ID == id)
      cnt <- called_counts[rws]

      if (all(cnt == 0)) {
        keep[rws] <- FALSE
        dropped_empty_ids <- c(dropped_empty_ids, id)
        emit("  %-14s DROPPED entirely -- all %d rows empty (rows %s).",
             id, length(rws), paste(rws, collapse = ", "))
        next
      }

      winner <- rws[which.max(cnt)]
      keep[setdiff(rws, winner)] <- FALSE

      # Where two non-empty rows exist, say how far they actually agree: a high
      # conflict count means the re-run disagrees with the original, which is a
      # genotyping problem the call-rate rule does not resolve.
      others <- setdiff(rws, winner)
      detail <- ""
      nonempty_others <- others[called_counts[others] > 0]
      if (length(nonempty_others) > 0) {
        conflicts <- 0L
        shared <- 0L
        for (o in nonempty_others) {
          a <- unlist(conformed_df[winner, marker_cols_dwr])
          b <- unlist(conformed_df[o, marker_cols_dwr])
          both <- !is.na(a) & !is.na(b)
          shared <- shared + sum(both)
          conflicts <- conflicts + sum(both & a != b)
        }
        detail <- sprintf(" | vs other row(s): %d shared calls, %d conflicting",
                          shared, conflicts)
        if (shared > 0 && conflicts > max_replicate_discordance * shared) {
          note("Duplicate %s: re-runs disagree at %d of %d shared calls (%.1f%%).",
               id, conflicts, shared, 100 * conflicts / shared)
        }
      }
      emit("  %-14s kept row %d (%d calls); dropped row(s) %s (%s calls)%s",
           id, winner, cnt[which.max(cnt)],
           paste(others, collapse = ", "),
           paste(called_counts[others], collapse = ", "), detail)
    }

    conformed_df <- conformed_df[keep, , drop = FALSE]
    rownames(conformed_df) <- NULL
    emit("[INFO] %d rows removed; %d rows remain.", sum(!keep), nrow(conformed_df))
    if (length(dropped_empty_ids) > 0) {
      emit("[WARNING] %d SAMPLE_IDs dropped entirely as every row was empty:",
           length(dropped_empty_ids))
      emit_obj(dropped_empty_ids)
    }
  } else if (any(duplicated(conformed_df$SAMPLE_ID))) {
    note("%d duplicate SAMPLE_IDs left in place (resolve_duplicate_ids = FALSE).",
         sum(duplicated(conformed_df$SAMPLE_ID)))
  }

  # --- Step 8: QC report -----------------------------------------------------
  rule <- strrep("=", 42)
  emit("")
  emit(rule)
  emit("              QC REPORT                   ")
  emit(rule)

  if (identical(colnames(conformed_df), dwr_cols)) {
    emit("[OK] 100% column alignment: output matches DWR reference structure exactly.")
  } else {
    emit("[ERROR] Column mismatch against DWR reference!")
  }

  emit("[INFO] Samples: %d in, %d out", n_samples, nrow(conformed_df))
  emit("[INFO] Markers matched: %d / %d", length(matches), length(marker_cols_dwr))

  missing_in_cfs <- setdiff(marker_cols_dwr, matches)
  if (length(missing_in_cfs) > 0) {
    emit("[INFO] %d DWR columns absent from CFS, filled with NA (e.g. Greb1l).",
         length(missing_in_cfs))
  }

  call_rates <- rowMeans(!is.na(conformed_df[, marker_cols_dwr, drop = FALSE]))
  names(call_rates) <- conformed_df$SAMPLE_ID
  emit("[INFO] Average sample call rate: %.2f%%", mean(call_rates) * 100)
  if (sum(call_rates < low_call_rate) > 0) {
    emit("[WARNING] %d samples have a call rate below %.0f%%.",
         sum(call_rates < low_call_rate), low_call_rate * 100)
  }
  if (sum(call_rates == 0) > 0) {
    emit("[WARNING] %d samples have no calls at all.", sum(call_rates == 0))
  }
  if (length(offset_samples) > 0) {
    emit("[WARNING] %d samples carried a shifted block (repair %s).",
         length(offset_samples), if (repair_offsets) "APPLIED" else "NOT applied")
  }
  if (length(dropped_empty_ids) > 0) {
    emit("[WARNING] %d SAMPLE_IDs dropped as every row was empty.",
         length(dropped_empty_ids))
  }

  if (any(duplicated(conformed_df$SAMPLE_ID))) {
    dups <- unique(conformed_df$SAMPLE_ID[duplicated(conformed_df$SAMPLE_ID)])
    emit("[WARNING] %d duplicate SAMPLE_IDs remain: %s",
         length(dups), paste(utils::head(dups, 5), collapse = ", "))
  } else {
    emit("[OK] SAMPLE_IDs are unique.")
  }

  if (length(gutted_loci) > 0) {
    emit("")
    emit("--- LOCI EMPTIED BY QC ---")
    emit("%d loci lost more than half their calls:", length(gutted_loci))
    emit_obj(gutted_loci)
    if (!repair_offsets && length(offset_samples) > 0) {
      emit("Most of this is the unrepaired offset above. Setting repair_offsets = TRUE")
      emit("recovers these rather than discarding them.")
    }
  }

  # Heterozygote excess: catches paralogous loci collapsed into one call, which
  # look perfectly valid to the width and allele checks.
  het <- paralog_scan(conformed_df, min_n = min_het_samples,
                      threshold = max_het_excess)
  ht <- het[het$flagged, ]
  if (nrow(ht) > 0) {
    emit("")
    emit("--- HETEROZYGOTE EXCESS ---")
    emit("%d loci are heterozygous far more often than their own allele",
         nrow(ht))
    emit("frequencies allow. This is what collapsed paralogs look like; the calls")
    emit("are individually well-formed, so nothing else here will flag them.")
    emit_obj(format(as.data.frame(
      ht[, c("locus", "n", "obs_het", "exp_het", "excess", "alleles")]
    ), digits = 3), row.names = FALSE)
    note("%d loci show heterozygote excess above %.2f (possible collapsed paralogs): %s",
         nrow(ht), max_het_excess, paste(ht$locus, collapse = ", "))
  }

  # Verify the CKMRsim invariant actually holds on the way out.
  bad_len_loci <- character(0)
  for (loc in unique(sub("\\.[12]$", "", marker_cols_dwr))) {
    c1 <- paste0(loc, ".1")
    c2 <- paste0(loc, ".2")
    if (!(c1 %in% dwr_cols && c2 %in% dwr_cols)) next
    a <- c(conformed_df[[c1]], conformed_df[[c2]])
    a <- a[.is_called(a)]
    if (length(unique(nchar(a))) > 1) bad_len_loci <- c(bad_len_loci, loc)
  }
  emit("")
  if (length(bad_len_loci) > 0) {
    emit("[WARNING] %d loci still have variable allele widths (CKMRsim will object):",
         length(bad_len_loci))
    emit_obj(utils::head(bad_len_loci, 10))
  } else {
    emit("[OK] All loci have uniform allele width (CKMRsim invariant satisfied).")
  }

  # What the run could NOT check. Unverified data must never read as verified.
  validated <- intersect(sub("\\.[12]$", "", matches), names(vcf_alleles))
  unverified <- sort(unique(setdiff(sub("\\.[12]$", "", matches), names(vcf_alleles))))
  emit("")
  emit("--- COVERAGE ---")
  emit("Loci validated against the VCF: %d", length(validated))
  if (length(unverified) > 0) {
    emit("Loci with NO VCF entry (width checked against the DWR CSV only, alleles unchecked): %d",
         length(unverified))
    emit_obj(utils::head(unverified, 20))
  }
  if (length(missing_in_cfs) > 0) {
    emit("DWR columns with no CFS source (empty in output): %d", length(missing_in_cfs))
  }

  emit("")
  emit("--- ANOMALIES ---")
  if (length(anomalies) == 0) {
    emit("None. Every check the run knows how to make came back clean.")
    emit("That is not proof the data is sound -- only that no known fault fired.")
  } else {
    emit("%d condition(s) need a human decision:", length(anomalies))
    for (i in seq_along(anomalies)) emit("  %d. %s", i, anomalies[i])
  }
  emit(rule)
  emit("")

  report <- structure(
    list(
      text = lines,
      anomalies = anomalies,
      offsets = offset_table,
      het_excess = het,
      call_rates = call_rates,
      coverage = list(
        validated = unique(validated),
        unverified = unverified,
        missing_in_cfs = missing_in_cfs
      )
    ),
    class = "conform_report"
  )

  out <- tibble::as_tibble(conformed_df)
  structure(
    out,
    class = c("geno_conformed", class(out)),
    conform_report = report
  )
}

#' The QC report from a conversion run
#'
#' @param x A `geno_conformed` object from [conform_genotypes()].
#' @return A `conform_report`: a list of `text` (every line of the report, in
#'   order), `anomalies` (conditions needing a human decision), `offsets` (the
#'   per-locus shift summary, or `NULL`), `het_excess` (the [paralog_scan()]
#'   result), `call_rates` (per sample, named by ID), and `coverage` (which
#'   loci could and could not be checked against a VCF). Printing it re-renders
#'   the report.
#' @seealso [conform_genotypes()]
#' @examples
#' syn <- function(f) {
#'   system.file("extdata", "synthetic", f, package = "genoCompaR")
#' }
#' out <- conform_genotypes(
#'   syn("synthetic_cfs_offset.csv"),
#'   syn("synthetic_dwr_reference.csv"),
#'   syn(c("synthetic_full.vcf", "synthetic_panel.vcf")),
#'   verbose = FALSE
#' )
#' conform_report(out)$anomalies
#' head(conform_report(out)$offsets)
#'
#' @export
conform_report <- function(x) {
  report <- attr(x, "conform_report")
  if (is.null(report)) {
    stop("`x` must be a geno_conformed object (see conform_genotypes()).",
         call. = FALSE)
  }
  report
}

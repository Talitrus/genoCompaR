# Reference-VCF-driven allele definitions.
#
# The central QC idea of conform_genotypes(): a VCF gives, per locus, both the
# number of SNP sites (= haplotype width) and the base permitted at each site,
# so an allele can be checked for validity outright instead of only for
# plausible width. That matters because a CFS column offset is only sometimes
# visible as a width change - when a locus and its neighbour have the same
# number of sites, a length test sees nothing wrong.

# Default reporting hooks, so these functions are usable on their own.
# conform_genotypes() passes its own collector and emitter instead.
.warn_note <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)
.silent_emit <- function(fmt, ...) invisible(NULL)

# Read reference VCFs into, for each locus, the set of bases permitted at each
# SNP site in order. A haplotype is valid iff its length equals the number of
# sites and every base is permitted at its position.
#
# `note` reports conditions needing a human decision; `emit` reports progress.
.load_vcf_alleles <- function(paths, note = .warn_note, emit = .silent_emit) {
  recs <- list()
  for (p in paths) {
    if (!file.exists(p)) {
      # Not fatal - the remaining references still validate what they cover -
      # but it silently shrinks QC coverage, so it must be surfaced.
      note("VCF not found, skipping: %s (allele checks lose coverage)", p)
      next
    }
    lines <- readLines(p, warn = FALSE)
    lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
    if (length(lines) == 0) next
    f <- strsplit(lines, "\t", fixed = TRUE)
    recs[[length(recs) + 1]] <- data.frame(
      locus = vapply(f, `[`, character(1), 1L),
      pos = as.integer(vapply(f, `[`, character(1), 2L)),
      ref = vapply(f, `[`, character(1), 4L),
      alt = vapply(f, `[`, character(1), 5L),
      stringsAsFactors = FALSE
    )
  }
  if (length(recs) == 0) {
    stop(
      "No VCF records loaded from: ", paste(paths, collapse = ", "),
      call. = FALSE
    )
  }
  v <- do.call(rbind, recs)
  # The panel VCF repeats loci already present in the full VCF; keep one row
  # per site, or the width doubles and every real allele turns invalid.
  v <- v[!duplicated(v[, c("locus", "pos")]), ]
  v <- v[order(v$locus, v$pos), ]

  # An indel cannot contribute one character to a fixed-width haplotype, so
  # sites with a multi-base REF are dropped entirely rather than counted - they
  # would otherwise inflate the expected width and admit multi-base "alleles".
  n_indel <- sum(nchar(v$ref) != 1L)
  if (n_indel > 0) {
    emit(
      "[INFO] Ignoring %d non-SNP (indel) VCF site(s); they cannot form fixed-width haplotypes.",
      n_indel
    )
  }
  v <- v[nchar(v$ref) == 1L, , drop = FALSE]
  if (nrow(v) == 0) {
    stop("No single-base VCF sites found; cannot define haplotypes.", call. = FALSE)
  }

  lapply(split(v, v$locus), function(d) {
    lapply(seq_len(nrow(d)), function(k) {
      alt <- strsplit(d$alt[k], ",", fixed = TRUE)[[1]]
      # Multi-base alternates at a SNP site are indel alleles: not permitted
      # values in the haplotype string.
      unique(c(d$ref[k], alt[nchar(alt) == 1L]))
    })
  })
}

# TRUE where the haplotype is a permitted allele of `allow` (a list of
# per-site character vectors). Uncalled values are FALSE, never NA.
.valid_against <- function(v, allow) {
  ok <- .is_called(v) & nchar(v) == length(allow)
  if (!any(ok)) return(ok)
  for (k in seq_along(allow)) {
    idx <- which(ok)
    ok[idx] <- substr(v[idx], k, k) %in% allow[[k]]
  }
  ok
}

# Widen each locus's permitted alleles with those actually seen in the
# reference.
#
# The VCFs are panel definitions called across a cohort, not a single fish -
# some sites carry three or four distinct alleles, which no diploid individual
# could produce. Their ALT lists are therefore a good allele catalogue, but not
# a guaranteed-complete one: a handful of reference alleles are absent from
# them, and those are a mix of genuinely rare variants and genotyping miscalls
# that cannot be told apart on frequency alone.
#
# Width still comes from the VCF, which is a property of the assay. Alleles are
# widened rather than hard-filtered so that a rare real variant is not deleted
# as corruption. The cost is that a miscall present in the reference is also
# admitted - see paralog_scan(), which catches that class instead.
.build_permitted <- function(vcf_alleles, dwr_ref, emit = .silent_emit) {
  loci <- names(vcf_alleles)
  extra_total <- 0L
  out <- vector("list", length(loci))
  names(out) <- loci

  for (locus in loci) {
    allow <- vcf_alleles[[locus]]
    observed <- character(0)
    for (col in paste0(locus, c(".1", ".2"))) {
      if (!col %in% colnames(dwr_ref)) next
      x <- dwr_ref[[col]]
      observed <- c(observed, x[.is_called(x)])
    }
    observed <- unique(observed)
    # Only same-width alleles: a different width is a defect, not an allele.
    observed <- observed[nchar(observed) == length(allow)]
    extra <- observed[!.valid_against(observed, allow)]
    extra_total <- extra_total + length(extra)
    out[[locus]] <- list(allow = allow, observed = observed)
  }
  if (extra_total > 0) {
    emit(
      "[INFO] %d reference alleles absent from the VCF ALT lists were admitted (cohort-called VCF).",
      extra_total
    )
  }
  out
}

# TRUE where the haplotype is permitted at `locus` by the VCF or seen in the
# reference.
.valid_allele <- function(v, permitted) {
  if (is.null(permitted)) return(rep(FALSE, length(v)))
  .valid_against(v, permitted$allow) | (.is_called(v) & v %in% permitted$observed)
}

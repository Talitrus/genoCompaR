# Fixtures for the CFS -> DWR conversion and paralog scan.
#
# Everything here runs against inst/extdata/synthetic, which contains no real
# fish: 5 files, ~50 kB, reproducing in miniature every fault the conversion is
# built to catch. See that directory's README.md for the feature -> test map,
# and data-raw/make_synthetic_fixtures.R for the (seeded) generator.

syn_path <- function(...) {
  system.file("extdata", "synthetic", ..., package = "genoCompaR")
}

SYN_CFS <- syn_path("synthetic_cfs_clean.csv")
SYN_OFFSET <- syn_path("synthetic_cfs_offset.csv")
SYN_DWR <- syn_path("synthetic_dwr_reference.csv")
SYN_VCFS <- syn_path(c("synthetic_full.vcf", "synthetic_panel.vcf"))

# conform_genotypes()'s defaults, named here so tests can reason about the CFS
# file's own columns without re-typing them.
SYN_ID_COL <- "indiv"
SYN_META <- c("indiv", "repunit", "collection", "sample_type",
              "rosa_genotypes_str")

.syn_cache <- new.env(parent = emptyenv())

cached <- function(key, value) {
  if (is.null(.syn_cache[[key]])) .syn_cache[[key]] <- value
  .syn_cache[[key]]
}

syn_vcf <- function() cached("vcf", .load_vcf_alleles(SYN_VCFS))
syn_dwr <- function() cached("dwr", .read_chr(SYN_DWR))
syn_cfs <- function(path = SYN_CFS) .read_chr(path)
syn_permitted <- function() cached("perm", .build_permitted(syn_vcf(), syn_dwr()))

# A conversion run. The two default-configuration runs are shared, since most
# tests want one of them; anything passing a non-default setting gets a fresh
# run, because a cached result would have been taken under the old settings.
syn_run <- function(cfs = SYN_CFS, ...) {
  args <- list(...)
  if (length(args) == 0) {
    return(cached(
      basename(cfs),
      conform_genotypes(cfs, SYN_DWR, SYN_VCFS, verbose = FALSE)
    ))
  }
  do.call(conform_genotypes,
          c(list(cfs, SYN_DWR, SYN_VCFS, verbose = FALSE), args))
}

syn_log <- function(cfs = SYN_CFS, ...) {
  paste(conform_report(syn_run(cfs, ...))$text, collapse = "\n")
}

syn_anomalies <- function(cfs = SYN_CFS, ...) {
  conform_report(syn_run(cfs, ...))$anomalies
}

# Alleles that are not permitted at the locus they sit in. The conversion's
# central claim is that this is zero on the way out.
count_invalid_alleles <- function(df, permitted = syn_permitted()) {
  bad <- 0L
  for (locus in names(permitted)) {
    for (col in paste0(locus, c(".1", ".2"))) {
      if (!col %in% colnames(df)) next
      vals <- df[[col]]
      vals <- vals[.is_called(vals)]
      if (length(vals) == 0) next
      bad <- bad + sum(!.valid_allele(vals, permitted[[locus]]))
    }
  }
  bad
}

n_called <- function(df) sum(!is.na(df[, setdiff(colnames(df), "SAMPLE_ID")]))

# Loci whose alleles are not all one width - what CKMRsim objects to.
variable_width_loci <- function(df) {
  offenders <- character(0)
  for (locus in .split_genotype_columns(colnames(df))$loci) {
    a <- unlist(df[, paste0(locus, c(".1", ".2"))], use.names = FALSE)
    a <- a[.is_called(a)]
    if (length(a) > 0 && length(unique(nchar(a))) > 1) {
      offenders <- c(offenders, locus)
    }
  }
  offenders
}

# Write a modified copy of a fixture to a temp file, for the failure-mode
# tests. `envir` keeps the file alive for the calling test only.
temp_csv <- function(df, envir = parent.frame()) {
  path <- withr::local_tempfile(fileext = ".csv", .local_envir = envir)
  utils::write.csv(df, path, row.names = FALSE, na = "")
  path
}

#!/usr/bin/env Rscript
# =============================================================================
# Generate the synthetic fixture set for the CFS -> DWR conversion.
#
#   Rscript data-raw/make_synthetic_fixtures.R
#
# Everything the test suite currently asserts against the real WY26 / FRH2024 /
# VCF files is reproduced here in miniature, on invented loci and invented fish.
# The real files are ~13 MB of genotypes belonging to CDFW/DWR; these are ~50 kB
# of nothing, so they can live in a package as inst/extdata and be shipped,
# regenerated and diffed freely.
#
# Files written to inst/extdata/synthetic/:
#
#   synthetic_full.vcf            SNP panel, plus indel records and loci that
#                                 the panel VCF repeats
#   synthetic_panel.vcf           second VCF, overlapping the first
#   synthetic_dwr_reference.csv   DWR-format reference genotypes
#   synthetic_cfs_clean.csv       CFS export, correct header order
#   synthetic_cfs_offset.csv      the same data with the 20260313 header fault
#
# The generator is seeded, so the outputs are stable across runs. Changing
# anything below changes the fixtures; tests/testthat/test-synthetic-fixtures.R
# is what proves the regenerated set still carries every property.
#
# See inst/extdata/synthetic/README.md for the feature -> test map.
# =============================================================================

set.seed(20260807)

this_file <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) sub("^--file=", "", hit[1]) else "data-raw/make_synthetic_fixtures.R"
})
ROOT <- normalizePath(file.path(dirname(this_file), ".."))
OUT  <- file.path(ROOT, "inst", "extdata", "synthetic")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

N_DWR_SAMPLES <- 80L
N_CFS_NORMAL  <- 32L
DWR_MISS_RATE <- 0.03
CFS_MISS_RATE <- 0.05

# -----------------------------------------------------------------------------
# Panel definition
#
# One entry per locus. `sites` is the VCF truth: a SNP site is a REF base plus
# ALT bases, and the haplotype width is the number of sites. `cfs` is the name
# the locus carries in a CFS export, `dwr` the name it carries in the DWR
# reference -- the pair is what transform_name() has to get right.
#
# Widths differ at every adjacent pair in the tag_id block, because that block
# is where the header-order fault is injected: a one-locus slip is only
# detectable where the neighbouring locus disagrees about what a valid genotype
# looks like.
# -----------------------------------------------------------------------------

snp <- function(ref, ...) list(ref = ref, alt = c(...))

# Bases the VCF permits at a site. Multi-base ALTs are indel alleles and cannot
# occupy one character of a fixed-width haplotype, so the loader drops them --
# tag_id_104 carries one to prove it.
permitted_bases <- function(site) unique(c(site$ref, site$alt[nchar(site$alt) == 1L]))

L <- function(dwr, cfs = NA_character_, sites = NULL, vcf = character(0),
              pool = NULL, pos = NULL) {
  list(dwr = dwr, cfs = cfs, sites = sites, vcf = vcf, pool = pool, pos = pos)
}

PANEL <- list(

  # --- NC_ block -------------------------------------------------------------
  # CFS writes genomic coordinates with dots throughout; DWR uses colon/dash.
  # The two flagged "both" appear in BOTH VCFs at identical (CHROM, POS): if the
  # loader does not dedupe, their widths double and every real allele at them
  # turns invalid.
  L("NC_037099.1:62937268-62937373", "NC_037099.1.62937268.62937373",
    sites = list(snp("A", "G"), snp("C", "T")),
    vcf = c("full", "panel"), pos = c(62L, 88L)),

  L("NC_037104.1:55923357-55923657", "NC_037104.1.55923357.55923657",
    sites = list(snp("A", "C"), snp("G", "T"), snp("C", "A")),
    vcf = "panel", pos = c(101L, 151L, 204L)),

  # Two indel records sit inside this locus in synthetic_full.vcf (added below).
  # Counted as sites they would make it 4 wide instead of 2.
  L("NC_037124.1:12281207-12281551", "NC_037124.1.12281207.12281551",
    sites = list(snp("G", "A"), snp("T", "C")),
    vcf = "full", pos = c(96L, 185L)),

  L("NC_037130.1:1062935-1063235", "NC_037130.1.1062935.1063235",
    sites = list(snp("T", "A"), snp("G", "A")),
    vcf = c("full", "panel"), pos = c(118L, 151L)),

  # --- Ots_ / OkiOts_ block --------------------------------------------------
  L("OkiOts_120255-113", "OkiOts_120255.113",
    sites = list(snp("A", "G"), snp("T", "C")),
    vcf = "full", pos = c(31L, 77L)),

  L("Ots_111312-435", "Ots_111312.435",
    sites = list(snp("C", "T")),
    vcf = "full", pos = 44L),

  # Every dot becomes a dash. CFS_to_DWR2.R replaced only the first one and
  # emitted "Ots_BMP-2.SNP1", which matches nothing in the DWR panel, so the
  # locus came out silently blank.
  L("Ots_BMP-2-SNP1", "Ots_BMP.2.SNP1",
    sites = list(snp("A", "G")),
    vcf = "full", pos = 60L),

  # The u07 family is the exception: only the FIRST dot becomes a dash, so the
  # DWR name keeps a literal dot in it.
  L("Ots_u07-07.161", "Ots_u07.07.161",
    sites = list(snp("G", "A"), snp("C", "T"), snp("A", "T")),
    vcf = "full", pos = c(22L, 51L, 109L)),

  L("Ots_u07-49.290", "Ots_u07.49.290",
    sites = list(snp("T", "C")),
    vcf = "full", pos = 35L),

  # --- tag_id_ block ---------------------------------------------------------
  # Identical names on both sides; only the allele suffix moves (CFS "x"/"x.1"
  # -> DWR "x.1"/"x.2"). This block carries the injected header fault.
  L("tag_id_101", "tag_id_101",
    sites = list(snp("A", "G"), snp("C", "T")),
    vcf = "full", pos = c(29L, 113L)),

  L("tag_id_102", "tag_id_102",
    sites = list(snp("G", "A"), snp("T", "C"), snp("A", "C")),
    vcf = "full", pos = c(10L, 72L, 94L)),

  L("tag_id_103", "tag_id_103",
    sites = list(snp("C", "T")),
    vcf = "full", pos = 63L),

  # Site 2 offers a multi-base ALT ("AT"). It must not become a permitted base.
  L("tag_id_104", "tag_id_104",
    sites = list(snp("A", "G"), snp("A", "G", "AT"), snp("T", "C"), snp("G", "A")),
    vcf = "full", pos = c(7L, 25L, 79L, 86L)),

  L("tag_id_105", "tag_id_105",
    sites = list(snp("T", "G"), snp("C", "A")),
    vcf = "full", pos = c(36L, 83L)),

  # The widest locus, and the one that carries the allele deliberately absent
  # from the VCF ALT lists (see NOVEL_ALLELE below).
  L("tag_id_106", "tag_id_106",
    sites = list(snp("A", "G"), snp("C", "T"), snp("A", "T"),
                 snp("G", "C"), snp("A", "C")),
    vcf = "full", pos = c(12L, 40L, 61L, 98L, 140L)),

  L("tag_id_107", "tag_id_107",
    sites = list(snp("G", "C"), snp("A", "T")),
    vcf = "full", pos = c(18L, 122L)),

  L("tag_id_108", "tag_id_108",
    sites = list(snp("C", "A"), snp("T", "G"), snp("G", "T")),
    vcf = "full", pos = c(45L, 90L, 131L)),

  L("tag_id_109", "tag_id_109",
    sites = list(snp("A", "T")),
    vcf = "full", pos = 57L),

  # Same width as tag_id_101 (2), so a shift landing here cannot be caught by
  # width alone -- its permitted bases at site 1 are disjoint from tag_id_101's
  # instead. This is the case a length-only QC would miss.
  L("tag_id_110", "tag_id_110",
    sites = list(snp("C", "T"), snp("A", "G")),
    vcf = "full", pos = c(21L, 108L)),

  # --- a locus with no VCF entry at all --------------------------------------
  # Present in both CSVs, described by neither VCF: its width can only be taken
  # from the DWR reference and its alleles cannot be checked. It is what the
  # COVERAGE section of the report exists to declare.
  # Its width (3) is shared with neither neighbour, so the offset scan cannot
  # mistake it for a shifted copy of one of them and quietly extend a run
  # through the one locus it is unable to check.
  L("tag_id_999", "tag_id_999", pool = c("ACG", "ATG", "GCG", "GTG")),

  # --- DWR-only markers ------------------------------------------------------
  # Winter-run Greb1l diagnostics are not on the CFS panel, so these columns
  # must appear in the output and be entirely empty.
  L("winter-greb1l-diagnostic_09_12269030_151", pool = c("A", "G")),
  L("winter-greb1l-diagnostic_10_12269914_151", pool = c("C", "T"))
)
names(PANEL) <- vapply(PANEL, `[[`, character(1), "dwr")

# A real DWR allele that no VCF ALT list mentions. The VCFs carry one sample
# column, so their ALT lists record what one fish happened to carry, not the
# population's alleles. build_permitted() has to widen the permitted set with
# what the DWR reference actually shows, or genuine variation reads as
# corruption -- and NULL_INVALID_ALLELES would then delete it.
NOVEL_LOCUS  <- "tag_id_106"
NOVEL_ALLELE <- "TCAGA"          # site 1 is 'T'; the VCF permits only A / G

# Loci that are in the VCFs, in panel order.
vcf_loci <- Filter(function(x) length(x$vcf) > 0, PANEL)

stopifnot(nchar(NOVEL_ALLELE) == length(PANEL[[NOVEL_LOCUS]]$sites))


# -----------------------------------------------------------------------------
# Allele pools
# -----------------------------------------------------------------------------

# Every haplotype the VCF permits at a locus, in site order.
all_haplotypes <- function(sites) {
  out <- ""
  for (s in sites) out <- as.vector(outer(out, permitted_bases(s), paste0))
  out
}

# Real loci segregate a handful of haplotypes, not the whole cross-product.
pool_for <- function(entry, n_max = 6L) {
  if (!is.null(entry$pool)) return(entry$pool)
  all <- sort(all_haplotypes(entry$sites))
  if (length(all) <= n_max) all else sort(sample(all, n_max))
}

POOLS <- lapply(PANEL, pool_for)
POOLS[[NOVEL_LOCUS]] <- c(POOLS[[NOVEL_LOCUS]], NOVEL_ALLELE)

# Skewed frequencies, so genotypes look like genotypes rather than uniform noise.
FREQS <- lapply(POOLS, function(p) {
  f <- sort(runif(length(p), 0.05, 1), decreasing = TRUE)
  f / sum(f)
})
# The novel allele is placed by hand, not drawn: it has to appear often enough
# to be visible and rarely enough to stay under 1% of all reference alleles.
FREQS[[NOVEL_LOCUS]] <- local({
  f <- FREQS[[NOVEL_LOCUS]]
  f[length(f)] <- 0
  f / sum(f)
})

# A genotype table: n rows, two columns per locus, in the given locus order.
draw_table <- function(n, loci, miss_rate) {
  m <- matrix(NA_character_, nrow = n, ncol = 2L * length(loci))
  for (k in seq_along(loci)) {
    pool <- POOLS[[loci[k]]]
    freq <- FREQS[[loci[k]]]
    a1 <- sample(pool, n, replace = TRUE, prob = freq)
    a2 <- sample(pool, n, replace = TRUE, prob = freq)
    drop <- runif(n) < miss_rate
    a1[drop] <- NA_character_
    a2[drop] <- NA_character_
    m[, 2L * k - 1L] <- a1
    m[, 2L * k]      <- a2
  }
  m
}

# Blank a random subset of loci in one row, down to roughly `keep` call rate.
thin_row <- function(row, keep) {
  n_loci <- length(row) %/% 2L
  drop <- sample(n_loci, round((1 - keep) * n_loci))
  for (k in drop) row[c(2L * k - 1L, 2L * k)] <- NA_character_
  row
}


# -----------------------------------------------------------------------------
# 1. VCFs
# -----------------------------------------------------------------------------

VCF_HEADER <- c(
  "##fileformat=VCFv4.2",
  "##source=synthetic fixture for CFS_to_DWR_haplotype; not real data",
  "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
  "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total alleles\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
          "FORMAT", "chinook2911"), collapse = "\t")
)

vcf_line <- function(chrom, pos, ref, alt, qual = "4211.7", gt = "0/1", dp = 31L) {
  paste(chrom, pos, ".", ref, alt, qual, ".",
        sprintf("AC=1;AN=2;DP=%d", dp), "GT:DP",
        sprintf("%s:%d", gt, dp), sep = "\t")
}

snp_records <- function(which_vcf) {
  out <- character(0)
  for (e in vcf_loci) {
    if (!(which_vcf %in% e$vcf)) next
    for (k in seq_along(e$sites)) {
      s <- e$sites[[k]]
      out <- c(out, vcf_line(e$dwr, e$pos[k], s$ref, paste(s$alt, collapse = ",")))
    }
  }
  out
}

# Indel records. Two sit inside a locus that also has SNP sites, so counting
# them would inflate that locus from 2 wide to 4 and admit multi-base
# "alleles"; the other two belong to a locus that is nothing BUT indels, and so
# has to vanish from the panel entirely rather than appear zero-wide.
# Lower-case bases mirror how the real full.vcf writes its Greb1l-region indels.
INDEL_RECORDS <- c(
  vcf_line("NC_037124.1:12281207-12281551", 213L, "taaaaaaaa", "taaaaaaa"),
  vcf_line("NC_037124.1:12281207-12281551", 241L, "AC", "A"),
  vcf_line("NC_037124.1:12270118-12270418", 119L,
           "tacacacacacacacacacacac", "tacacacacacacacacac,tacacacacacacacac"),
  vcf_line("NC_037124.1:12270118-12270418", 158L, "GTT", "G")
)

writeLines(c(VCF_HEADER, snp_records("full"), INDEL_RECORDS),
           file.path(OUT, "synthetic_full.vcf"))
writeLines(c(VCF_HEADER, snp_records("panel")),
           file.path(OUT, "synthetic_panel.vcf"))


# -----------------------------------------------------------------------------
# 2. DWR reference
#
# Column order is sorted by locus name, which is NOT the order the CFS export
# uses -- the conversion has to reorder, not just rename.
# -----------------------------------------------------------------------------

dwr_loci <- sort(vapply(PANEL, `[[`, character(1), "dwr"))
dwr_cols <- c("SAMPLE_ID", as.vector(rbind(paste0(dwr_loci, ".1"),
                                           paste0(dwr_loci, ".2"))))

dwr_geno <- draw_table(N_DWR_SAMPLES, dwr_loci, DWR_MISS_RATE)

# Seed the novel allele into a few reference fish. It has to be present, or
# the widening step has nothing to widen with and the fixture stops testing it.
novel_col <- 2L * match(NOVEL_LOCUS, dwr_loci) - 1L
dwr_geno[c(4L, 17L, 33L, 58L, 71L), novel_col] <- NOVEL_ALLELE

dwr <- data.frame(SAMPLE_ID = sprintf("FRHSYN%03d", seq_len(N_DWR_SAMPLES)),
                  dwr_geno, stringsAsFactors = FALSE)
colnames(dwr) <- dwr_cols

write.csv(dwr, file.path(OUT, "synthetic_dwr_reference.csv"),
          row.names = FALSE, quote = FALSE, na = "NA")


# -----------------------------------------------------------------------------
# 3. CFS export (clean)
#
# Column order here is the CFS one: metadata block, then NC_, Ots_/OkiOts_ and
# tag_id_ marker pairs in the order PANEL declares them.
# -----------------------------------------------------------------------------

cfs_entries <- Filter(function(x) !is.na(x$cfs), PANEL)
cfs_loci    <- vapply(cfs_entries, `[[`, character(1), "dwr")   # PANEL order
cfs_names   <- vapply(cfs_entries, `[[`, character(1), "cfs")
cfs_marker_cols <- as.vector(rbind(cfs_names, paste0(cfs_names, ".1")))

# --- the sample roster -------------------------------------------------------
# Every row is one of: an ordinary fish, a deliberate call-rate case, or one of
# the four duplicate-ID shapes the resolution rule has to handle.
roster <- rbind(
  data.frame(id = sprintf("SYN%04d", seq_len(N_CFS_NORMAL)),
             role = "normal", keep = 1, stringsAsFactors = FALSE),

  # Call-rate reporting: one poor sample, one that failed outright. Both have
  # unique IDs, so neither is touched by duplicate resolution.
  data.frame(id = "SYN_LOWCALL", role = "thin",   keep = 0.35),
  data.frame(id = "SYN_NOCALL",  role = "empty",  keep = 0),

  # The ordinary re-run: a failed first attempt and a good second one. Keep the
  # row with more calls.
  data.frame(id = "SYN_DUP_HIGHLOW", role = "dup_low",  keep = 0.30),
  data.frame(id = "SYN_DUP_HIGHLOW", role = "dup_high", keep = 0.95),

  # Both attempts failed. The ID contributes nothing and is dropped entirely.
  data.frame(id = "SYN_DUP_EMPTY", role = "empty", keep = 0),
  data.frame(id = "SYN_DUP_EMPTY", role = "empty", keep = 0),

  # Two good runs of one fish that disagree. Call rate picks a row; it does not
  # explain the disagreement, so this must surface as an anomaly.
  data.frame(id = "SYN_DUP_DISCORD", role = "discord_a", keep = 0.88),
  data.frame(id = "SYN_DUP_DISCORD", role = "discord_b", keep = 0.90),

  # Three rows for one ID: one dead, two live. Only the best survives.
  data.frame(id = "SYN_DUP_TRIPLE", role = "empty",    keep = 0),
  data.frame(id = "SYN_DUP_TRIPLE", role = "dup_low",  keep = 0.60),
  data.frame(id = "SYN_DUP_TRIPLE", role = "dup_high", keep = 0.90)
)
n_cfs <- nrow(roster)

cfs_geno <- draw_table(n_cfs, cfs_loci, CFS_MISS_RATE)

# Rows sharing a SAMPLE_ID are re-runs of ONE fish, so they start from one
# genotype and differ only in how much of it the run managed to call. Drawing
# them independently would make every duplicate look like a replicate conflict
# and bury the one pair that genuinely is one.
for (id in unique(roster$id[duplicated(roster$id)])) {
  rws <- which(roster$id == id)
  for (r in rws[-1]) cfs_geno[r, ] <- cfs_geno[rws[1], ]
}

# The exception: two good runs of one fish that disagree at ~15% of loci. The
# disagreeing values are still legitimate alleles of their locus -- a replicate
# conflict is not a corrupt allele, and the QC must not confuse the two.
ib <- which(roster$role == "discord_b")
for (k in sample(length(cfs_loci), round(0.30 * length(cfs_loci)))) {
  pool <- POOLS[[cfs_loci[k]]]
  cfs_geno[ib, 2L * k - 1L] <- sample(setdiff(pool, cfs_geno[ib, 2L * k - 1L]), 1)
}

# Call rates last, so a thinned row is a partial read of its own genotype.
for (i in seq_len(n_cfs)) {
  if (roster$keep[i] == 0) {
    cfs_geno[i, ] <- NA_character_
  } else if (roster$keep[i] < 1) {
    cfs_geno[i, ] <- thin_row(cfs_geno[i, ], roster$keep[i])
  }
}

# The novel allele turns up in incoming data too, not only in the reference.
# Without widening these would be flagged as unrecognised.
novel_cfs_col <- 2L * match(NOVEL_LOCUS, cfs_loci) - 1L
cfs_geno[c(3L, 21L), novel_cfs_col] <- NOVEL_ALLELE

# A couple of genuinely empty cells rather than the literal string NA, since
# exports have carried both and is_called() has to treat them alike.
cfs_geno[9L, 3:4] <- ""

cfs <- data.frame(
  indiv       = roster$id,
  repunit     = "NA",
  collection  = ifelse(grepl("^dup_(high|low)$|^discord_b$", roster$role),
                       "20260315_CVP_Rapid_Order2044", "Mhap_Order2001_i11i12"),
  sample_type = "mixed",
  rosa_genotypes_str = vapply(seq_len(n_cfs), function(i)
    paste(sample(c("E", "N"), 12, TRUE), collapse = ""), character(1)),
  cfs_geno,
  stringsAsFactors = FALSE
)
colnames(cfs) <- c("indiv", "repunit", "collection", "sample_type",
                   "rosa_genotypes_str", cfs_marker_cols)

write.csv(cfs, file.path(OUT, "synthetic_cfs_clean.csv"),
          row.names = FALSE, quote = FALSE, na = "NA")


# -----------------------------------------------------------------------------
# 4. CFS export (header-order fault)
#
# This reproduces the 20260313 export exactly as it failed: the genotype
# columns come off the instrument in the right order, but the header lists the
# first locus of the tag_id block at the END of the block. Every column after
# it is therefore labelled with its neighbour's name, and the mislabelling is
# invisible to anything that only checks widths -- which is the whole reason
# the offset scan validates alleles against the VCF instead.
#
# Data is untouched; only the header line is permuted. Reading the result: the
# value under label `tag_id_102` is really tag_id_101's genotype, i.e. it
# belongs one locus to the LEFT, which the scan reports as shift -1.
# -----------------------------------------------------------------------------

hdr <- colnames(cfs)
blk <- match(c("tag_id_101", "tag_id_101.1"), hdr)       # the pair that moved
end <- match("tag_id_110.1", hdr)                        # end of the tag block
stopifnot(!anyNA(c(blk, end)), blk[2] == blk[1] + 1L, end > blk[2])

moved <- append(hdr[-blk], hdr[blk], after = end - 2L)
stopifnot(setequal(moved, hdr), !identical(moved, hdr), length(moved) == length(hdr))

broken <- cfs
colnames(broken) <- moved
write.csv(broken, file.path(OUT, "synthetic_cfs_offset.csv"),
          row.names = FALSE, quote = FALSE, na = "NA")


# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

n_sites <- sum(vapply(vcf_loci, function(e) length(e$sites), integer(1)))
cat(sprintf("Wrote 5 files to %s\n", OUT))
cat(sprintf("  loci: %d total | %d in the VCFs | %d on the CFS panel | %d DWR-only\n",
            length(PANEL), length(vcf_loci), length(cfs_entries),
            sum(is.na(vapply(PANEL, `[[`, character(1), "cfs")))))
cat(sprintf("  VCF SNP sites: %d (+%d indel records that must be ignored)\n",
            n_sites, length(INDEL_RECORDS)))
cat(sprintf("  DWR reference: %d samples x %d columns\n", nrow(dwr), ncol(dwr)))
cat(sprintf("  CFS export:    %d rows (%d distinct IDs) x %d columns\n",
            nrow(cfs), length(unique(cfs$indiv)), ncol(cfs)))
cat(sprintf("  haplotype widths present: %s\n",
            paste(sort(unique(vapply(vcf_loci, function(e) length(e$sites),
                                     integer(1)))), collapse = ", ")))

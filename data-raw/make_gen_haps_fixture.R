#!/usr/bin/env Rscript
# =============================================================================
# Generate the synthetic fixture for gen_haps().
#
#   Rscript data-raw/make_gen_haps_fixture.R
#
# Writes inst/extdata/synthetic/synthetic_unfiltered_haplotypes.csv: a raw,
# long-format "*_observed_unfiltered_haplotype.csv" table (one row per allele
# call) small enough to ship, that deliberately exercises every stage of
# gen_haps() at the pipeline's default thresholds (haplotype_depth = 4,
# total_depth = 8, allele_balance = 0.35). Every row below is placed by hand,
# not sampled - each one targets exactly one stage of the pipeline, so the
# expected output can be worked out by inspection rather than by re-running a
# reference implementation. See inst/extdata/synthetic/README.md for the
# feature -> test map and the worked-out expected result.
#
# Two clean loci (locA, locB) carry the baseline het/hom calls and stay
# present in the final output; three fault loci (locC, locD, locE) each
# belong to exactly one sample and are entirely eliminated by their fault, so
# they vanish from the output altogether and exercise the "loci dropped"
# anomaly.
# =============================================================================

this_file <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) sub("^--file=", "", hit[1]) else "data-raw/make_gen_haps_fixture.R"
})
ROOT <- normalizePath(file.path(dirname(this_file), ".."))
OUT  <- file.path(ROOT, "inst", "extdata", "synthetic")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

row_ <- function(indiv, group, locus, haplo, depth, rank, ab = 1.0) {
  data.frame(
    group = group, indiv.ID = indiv, locus = locus, haplo = haplo,
    depth = depth, sum.Phred.C = depth - 0.5, max.Phred.C = 0.999,
    allele.balance = ab, rank = rank, stringsAsFactors = FALSE
  )
}

rows <- rbind(
  # --- FISH01: baseline, everything a clean het/hom pair -------------------
  row_("FISH01", "gA", "locA", "AC", 50, 1),
  row_("FISH01", "gA", "locA", "AG", 45, 2),
  row_("FISH01", "gA", "locB", "CC", 60, 1),   # homozygous single call

  # --- FISH02: baseline, a second normal sample for depth summaries --------
  row_("FISH02", "gA", "locA", "AC", 40, 1),
  row_("FISH02", "gA", "locA", "AA", 38, 2),
  row_("FISH02", "gA", "locB", "CC", 55, 1),
  row_("FISH02", "gA", "locB", "CT", 50, 2),

  # --- FISH03: allele_balance filter (0.20 < 0.35) --------------------------
  # locA's only call is dropped before NAlleles/TotalDepth are even computed,
  # so locA is NA for FISH03 while staying present for everyone else.
  row_("FISH03", "gA", "locA", "AC", 40, 1, ab = 0.20),
  row_("FISH03", "gA", "locB", "CC", 42, 1),

  # --- FISH04: haplotype_depth filter (2 < 4) -------------------------------
  row_("FISH04", "gA", "locA", "AC", 2, 1),
  row_("FISH04", "gA", "locB", "CC", 44, 1),

  # --- FISH05: total_depth filter (5 < 8) -----------------------------------
  # A single call whose OWN depth clears haplotype_depth but whose (here,
  # single-row) TotalDepth does not - the case per-call depth filtering alone
  # cannot catch.
  row_("FISH05", "gA", "locE", "GT", 5, 1),
  row_("FISH05", "gA", "locA", "AC", 41, 1),
  row_("FISH05", "gA", "locB", "CC", 46, 1),

  # --- FISH06: N allele wipes the WHOLE sample-locus pair -------------------
  # locD's rank-1 call is perfectly good; it is dropped anyway because rank-2
  # at the same locus is "NN" - verified against the real pipeline (see
  # CLAUDE.md).
  row_("FISH06", "gA", "locD", "AT", 48, 1),
  row_("FISH06", "gA", "locD", "NN", 44, 2),
  row_("FISH06", "gA", "locA", "AC", 43, 1),
  row_("FISH06", "gA", "locB", "CC", 47, 1),

  # --- FISH07: >2 surviving alleles ("contaminated") locus ------------------
  row_("FISH07", "gA", "locC", "AA", 30, 1),
  row_("FISH07", "gA", "locC", "AC", 28, 2),
  row_("FISH07", "gA", "locC", "AG", 26, 3),
  row_("FISH07", "gA", "locA", "AC", 39, 1),
  row_("FISH07", "gA", "locB", "CC", 41, 1),

  # --- FISH08: every call filtered out -> lands in missing_samples ---------
  row_("FISH08", "gA", "locA", "AC", 1, 1),
  row_("FISH08", "gA", "locB", "CC", 1, 1),

  # --- DUP01: re-run under two group tags -----------------------------------
  # "runA" calls only locA (1 locus); "runB" calls locA + locB (2 loci) and
  # wins duplicate resolution, so DUP01's final row carries runB's calls.
  row_("DUP01", "runA", "locA", "AC", 44, 1),
  row_("DUP01", "runB", "locA", "AC", 42, 1),
  row_("DUP01", "runB", "locB", "CC", 40, 1)
)

write.csv(rows, file.path(OUT, "synthetic_unfiltered_haplotypes.csv"),
          row.names = FALSE, quote = FALSE, na = "NA")

cat(sprintf("Wrote %d rows (%d distinct IDs, %d loci) to %s\n",
            nrow(rows), length(unique(rows$indiv.ID)), length(unique(rows$locus)),
            file.path(OUT, "synthetic_unfiltered_haplotypes.csv")))

# gen_haps() ports the Nextflow chinook_gt pipeline's GEN_HAPS process
# (modules/microhaplot.nf -> bin/gen_haps.R) into dependency-free base R. The
# fixture (inst/extdata/synthetic/synthetic_unfiltered_haplotypes.csv,
# data-raw/make_gen_haps_fixture.R) hand-places one row set per pipeline
# stage, so every expected value here was worked out by inspection - not
# copied from a run - then cross-checked against the actual output. The
# algorithm itself (including the two non-obvious quirks called out in
# ?gen_haps) was separately verified against the real, installed
# microhaplotopia@emc-edits package on real shipped data
# (inst/extdata/AFTC_ILC_P134_full_filtered_haplotypes.csv): zero cell
# differences, identical row/column order.

SYN_HAPS <- system.file("extdata", "synthetic", "synthetic_unfiltered_haplotypes.csv",
                        package = "genoCompaR")

syn_gen_haps <- function(...) gen_haps(SYN_HAPS, verbose = FALSE, ...)

test_that("thresholds default to the Nextflow pipeline's values", {
  expect_equal(formals(gen_haps)$haplotype_depth, 4)
  expect_equal(formals(gen_haps)$total_depth, 8)
  expect_equal(formals(gen_haps)$allele_balance, 0.35)
})

test_that("invalid thresholds are rejected the way main.nf rejects them", {
  expect_error(gen_haps(SYN_HAPS, allele_balance = -0.1), "between 0 and 1")
  expect_error(gen_haps(SYN_HAPS, allele_balance = 1.1), "between 0 and 1")
  expect_error(gen_haps(SYN_HAPS, haplotype_depth = -1), "non-negative")
  expect_error(gen_haps(SYN_HAPS, total_depth = -1), "non-negative")
})

test_that("result shape: 9 samples, 2 surviving loci, group/indiv.ID first", {
  out <- syn_gen_haps()
  expect_s3_class(out, "geno_haps")
  expect_equal(nrow(out), 9L)
  expect_equal(names(out)[1:2], c("group", "indiv.ID"))
  expect_setequal(get_loci(as_geno_data(out)), c("locA", "locB"))
})

test_that("group is 'ND' on every row, not just the re-added ones", {
  # Verified quirk: the real pipeline does this unconditionally. See
  # ?gen_haps Details.
  out <- syn_gen_haps()
  expect_true(all(out$group == "ND"))
})

test_that("preserve_group = TRUE keeps each sample's own group instead of 'ND'", {
  out <- syn_gen_haps(preserve_group = TRUE)

  # DUP01 takes whichever group won duplicate resolution (Step 6), not "gA".
  expect_equal(out$group[out$indiv.ID == "DUP01"], "runB")

  # A normal survivor keeps its own (only) group.
  expect_equal(out$group[out$indiv.ID == "FISH01"], "gA")

  # A sample re-added in Step 8 (no calls survived) still keeps the group its
  # first raw call carried, even though it contributes no genotype data.
  expect_equal(out$group[out$indiv.ID == "FISH08"], "gA")
})

test_that("preserve_group defaults to FALSE and is echoed in the report params", {
  out <- syn_gen_haps()
  expect_false(gen_haps_report(out)$params$preserve_group)
  expect_true(gen_haps_report(syn_gen_haps(preserve_group = TRUE))$params$preserve_group)
})

test_that("a normal heterozygous call pairs both alleles", {
  out <- syn_gen_haps()
  fish01 <- out[out$indiv.ID == "FISH01", ]
  expect_equal(fish01$locA.1, "AC")
  expect_equal(fish01$locA.2, "AG")
})

test_that("a homozygous single call is duplicated into both allele columns", {
  out <- syn_gen_haps()
  fish01 <- out[out$indiv.ID == "FISH01", ]
  expect_equal(fish01$locB.1, "CC")
  expect_equal(fish01$locB.2, "CC")
})

test_that("allele_balance filters a single call, leaving that locus NA", {
  out <- syn_gen_haps()
  fish03 <- out[out$indiv.ID == "FISH03", ]
  expect_true(is.na(fish03$locA.1) && is.na(fish03$locA.2))
  # locB is untouched - the filter is per-call, not per-sample.
  expect_equal(fish03$locB.1, "CC")
})

test_that("haplotype_depth filters a single call, leaving that locus NA", {
  out <- syn_gen_haps()
  fish04 <- out[out$indiv.ID == "FISH04", ]
  expect_true(is.na(fish04$locA.1) && is.na(fish04$locA.2))
  expect_equal(fish04$locB.1, "CC")
})

test_that("total_depth catches what per-call depth filtering cannot", {
  # locE's one call passes haplotype_depth (5 >= 4) but fails total_depth
  # (5 < 8); the locus disappears from the output entirely as a result (it's
  # FISH05's only call there).
  out <- syn_gen_haps()
  expect_false("locE.1" %in% names(out))
  fish05 <- out[out$indiv.ID == "FISH05", ]
  expect_equal(fish05$locA.1, "AC")
  expect_equal(fish05$locA.2, "AC")
})

test_that("an N/X call wipes the whole sample-locus pair, not just itself", {
  # locD's rank-1 call ("AT") is individually well-formed and would pass every
  # filter on its own; it is dropped anyway because rank-2 at the same locus
  # is "NN". Since FISH06 is locD's only sample, the locus vanishes entirely.
  out <- syn_gen_haps()
  expect_false("locD.1" %in% names(out))

  report <- gen_haps_report(out)
  expect_equal(nrow(report$nxa), 1L)
  expect_equal(report$nxa$locus, "locD")
  expect_equal(report$nxa$haplo, "NN")
})

test_that("more than 2 surviving alleles drops the whole locus as contaminated", {
  out <- syn_gen_haps()
  expect_false("locC.1" %in% names(out))

  report <- gen_haps_report(out)
  expect_equal(nrow(report$xtraalleles), 3L)
  expect_true(all(report$xtraalleles$locus == "locC"))
  expect_equal(report$xtraalleles_individuals$count[report$xtraalleles_individuals$indiv.ID == "FISH07"], 1L)
  expect_equal(report$xtraalleles_locus$count[report$xtraalleles_locus$locus == "locC"], 1L)
})

test_that("loci with no surviving calls are dropped and noted as an anomaly", {
  out <- syn_gen_haps()
  report <- gen_haps_report(out)
  expect_false(any(c("locC.1", "locD.1", "locE.1") %in% names(out)))
  expect_true(any(grepl("locC, locD, locE", report$anomalies)))
})

test_that("a sample with no surviving calls at all is re-added as an all-NA row", {
  out <- syn_gen_haps()
  fish08 <- out[out$indiv.ID == "FISH08", ]
  expect_equal(nrow(fish08), 1L)
  expect_true(all(is.na(fish08[, c("locA.1", "locA.2", "locB.1", "locB.2")])))
  expect_equal(fish08$group, "ND")

  report <- gen_haps_report(out)
  expect_equal(report$missing_samples$indiv.ID, "FISH08")
  expect_equal(report$missing_samples$group, "gA")  # original group, pre-ND override
})

test_that("a duplicate sample ID keeps the run with the most loci called", {
  # DUP01/runA calls only locA (1 locus); DUP01/runB calls locA + locB
  # (2 loci) and wins. Verified against the real
  # microhaplotopia::resolve_duplicate_samples(resolve = "drop") on a
  # constructed tie/non-tie case during development.
  out <- syn_gen_haps()
  dup <- out[out$indiv.ID == "DUP01", ]
  expect_equal(nrow(dup), 1L)
  expect_equal(dup$locA.1, "AC")
  expect_equal(dup$locB.1, "CC")

  report <- gen_haps_report(out)
  expect_setequal(report$duplicate_samples$group, c("runA", "runB"))
})

test_that("as_geno_data() accepts the result with its defaults unchanged", {
  out <- syn_gen_haps()
  geno <- as_geno_data(out)
  expect_s3_class(geno, "geno_data")
  expect_equal(nrow(geno), nrow(out))
  expect_setequal(get_loci(geno), c("locA", "locB"))
})

test_that("the print method summarizes samples, loci and anomaly count", {
  out <- syn_gen_haps()
  txt <- paste(capture.output(print(out)), collapse = "\n")
  expect_match(txt, "<geno_haps> 9 samples x 2 loci")
  expect_match(txt, "3 anomalies")
})

test_that("gen_haps_report() round-trips and rejects non-geno_haps input", {
  out <- syn_gen_haps()
  report <- gen_haps_report(out)
  expect_s3_class(report, "gen_haps_report")
  expect_type(report$text, "character")
  expect_true(length(report$text) > 0)
  expect_equal(report$params, list(haplotype_depth = 4, total_depth = 8,
                                   allele_balance = 0.35, preserve_group = FALSE))

  expect_error(gen_haps_report(data.frame(x = 1)), "geno_haps")
})

test_that("thresholds are adjustable and change the result", {
  # Loosen allele_balance enough to let FISH03's 0.20 call through.
  out <- syn_gen_haps(allele_balance = 0.1)
  fish03 <- out[out$indiv.ID == "FISH03", ]
  expect_equal(fish03$locA.1, "AC")
  expect_equal(fish03$locA.2, "AC")
})

test_that("an already-loaded data frame works the same as a file path", {
  df <- utils::read.csv(SYN_HAPS, colClasses = "character", check.names = FALSE)
  out_df <- gen_haps(df, verbose = FALSE)
  out_path <- syn_gen_haps()
  expect_equal(as.data.frame(out_df), as.data.frame(out_path))
})

test_that("missing required columns are a clear error", {
  df <- utils::read.csv(SYN_HAPS, colClasses = "character", check.names = FALSE)
  df$rank <- NULL
  expect_error(gen_haps(df, verbose = FALSE), "rank")
})

test_that("a nonexistent file path is a clear error", {
  expect_error(gen_haps("/no/such/file.csv", verbose = FALSE), "not found")
})

# How the conversion behaves on inputs it was not designed for. The point is
# that it fails loudly on things it cannot handle and degrades visibly on
# things it can, rather than producing a plausible-looking file either way.

test_that("a missing input file fails with a usable message", {
  missing <- syn_path("nope.csv")
  expect_error(conform_genotypes(missing, SYN_DWR, SYN_VCFS, verbose = FALSE),
               "Required input not found")
  # The message names the working directory, since a relative path is the
  # usual cause.
  err <- tryCatch(
    conform_genotypes(missing, SYN_DWR, SYN_VCFS, verbose = FALSE),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "working directory")
})

test_that("a CFS file without the ID column fails rather than guessing", {
  cfs <- syn_cfs()
  names(cfs)[names(cfs) == SYN_ID_COL] <- "not_the_id"

  expect_error(
    conform_genotypes(temp_csv(cfs), SYN_DWR, SYN_VCFS, verbose = FALSE),
    "no '.*' column to use as SAMPLE_ID"
  )
})

test_that("a DWR reference without SAMPLE_ID fails", {
  dwr <- syn_dwr()
  names(dwr)[names(dwr) == "SAMPLE_ID"] <- "something_else"

  expect_error(
    conform_genotypes(SYN_CFS, temp_csv(dwr), SYN_VCFS, verbose = FALSE),
    "no SAMPLE_ID column"
  )
})

test_that("a group/indiv.ID reference header is accepted and noted", {
  # Some DWR exports lead with "group,indiv.ID,..." instead of SAMPLE_ID.
  dwr <- syn_dwr()
  names(dwr)[names(dwr) == "SAMPLE_ID"] <- "indiv.ID"
  dwr <- cbind(group = "syn", dwr, stringsAsFactors = FALSE)

  out <- conform_genotypes(SYN_CFS, temp_csv(dwr), SYN_VCFS, verbose = FALSE)
  expect_identical(colnames(out), colnames(syn_dwr()))
  expect_true(any(grepl("indiv.ID renamed to SAMPLE_ID",
                        conform_report(out)$anomalies)))
})

test_that("an empty CFS file fails", {
  expect_error(
    conform_genotypes(temp_csv(syn_cfs()[0, ]), SYN_DWR, SYN_VCFS, verbose = FALSE),
    "no data rows"
  )
})

test_that("missing optional metadata is noted but not fatal", {
  cfs <- syn_cfs()
  cfs$repunit <- NULL
  cfs$sample_type <- NULL

  out <- conform_genotypes(temp_csv(cfs), SYN_DWR, SYN_VCFS, verbose = FALSE)
  expect_equal(ncol(out), ncol(syn_dwr()))
  expect_true(any(grepl("metadata columns absent", conform_report(out)$anomalies)))
})

test_that("unpaired allele columns are raised as an anomaly", {
  cfs <- syn_cfs()
  # Drop one half of an allele pair.
  cfs[[setdiff(colnames(cfs), SYN_META)[1]]] <- NULL

  out <- conform_genotypes(temp_csv(cfs), SYN_DWR, SYN_VCFS, verbose = FALSE)
  expect_true(any(grepl("no partner|pair evenly", conform_report(out)$anomalies)))
})

test_that("losing a VCF degrades coverage visibly instead of silently", {
  # Drop the full VCF, leaving only the 3-locus panel.
  out <- conform_genotypes(SYN_CFS, SYN_DWR, SYN_VCFS[2], verbose = FALSE)
  txt <- paste(conform_report(out)$text, collapse = "\n")

  # Still produces a file, but says plainly what it could not check.
  expect_identical(colnames(out), colnames(syn_dwr()))
  expect_match(txt, "NO VCF entry")
  expect_match(txt, "alleles unchecked")
  expect_gt(length(conform_report(out)$coverage$unverified), 1L)
})

test_that("a VCF path that does not exist is an anomaly, not a failure", {
  out <- conform_genotypes(SYN_CFS, SYN_DWR, c(SYN_VCFS, syn_path("no_such.vcf")),
                           verbose = FALSE)
  expect_true(any(grepl("VCF not found", conform_report(out)$anomalies)))
  expect_equal(count_invalid_alleles(out), 0L)
})

test_that("anomalies belong to their own run, not to the session", {
  # The collector lives in the run, so a run that found a fault cannot leak it
  # into the next one - and re-reading an old result still shows what it found.
  broken <- syn_run(SYN_OFFSET)
  clean <- conform_genotypes(SYN_CFS, SYN_DWR, SYN_VCFS, verbose = FALSE)

  expect_true(any(grepl("repair is OFF", conform_report(broken)$anomalies)))
  expect_false(any(grepl("repair is OFF", conform_report(clean)$anomalies)))
  # ... and the earlier result is unchanged by the later run.
  expect_true(any(grepl("repair is OFF", conform_report(broken)$anomalies)))
})

test_that("a clean run reports nothing beyond the deliberate discordance", {
  # SYN_DUP_DISCORD is the one condition the fixture is built to raise.
  # Anything else appearing here is new and wants looking at.
  unexpected <- grep("re-runs disagree", syn_anomalies(), invert = TRUE,
                     value = TRUE)
  expect_equal(unexpected, character(0),
               info = paste("unexpected anomalies:",
                            paste(unexpected, collapse = " | ")))
})

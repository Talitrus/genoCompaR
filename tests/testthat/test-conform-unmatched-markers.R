# CFS markers that harmonize to a name absent from the DWR panel: kept by
# default (keep_unmatched_markers = TRUE), droppable for the old behaviour.

test_that("a CFS-only marker is appended beyond the DWR template by default", {
  cfs <- syn_cfs()
  cfs$tag_id_extra <- "AC"
  cfs$tag_id_extra.1 <- "AC"

  out <- conform_genotypes(temp_csv(cfs), SYN_DWR, SYN_VCFS, verbose = FALSE)

  expect_identical(colnames(out),
                   c(colnames(syn_dwr()), "tag_id_extra.1", "tag_id_extra.2"))
  expect_true(all(out$tag_id_extra.1 == "AC"))
  expect_true(all(out$tag_id_extra.2 == "AC"))

  txt <- paste(conform_report(out)$text, collapse = "\n")
  expect_match(txt, "kept as extra columns")
  expect_match(txt, "tag_id_extra")

  coverage <- conform_report(out)$coverage
  expect_true("tag_id_extra.1" %in% coverage$extra_cols)
})

test_that("a kept extra marker converts cleanly into the rest of the package", {
  cfs <- syn_cfs()
  cfs$tag_id_extra <- "AC"
  cfs$tag_id_extra.1 <- "AC"

  out <- conform_genotypes(temp_csv(cfs), SYN_DWR, SYN_VCFS, verbose = FALSE)
  geno <- as_geno_data(out, id_col = "SAMPLE_ID", group_col = NULL)

  expect_true("tag_id_extra" %in% get_loci(geno))
  expect_equal(length(get_loci(geno)), 23L)  # 22 DWR loci + 1 kept extra

  cmp <- compare_genotypes(geno, "SYN0001", "SYN0002", output = "tibble")
  expect_equal(cmp$n_loci_total[[1]], 23L)
})

test_that("keep_unmatched_markers = FALSE reproduces the old drop behaviour", {
  cfs <- syn_cfs()
  cfs$tag_id_extra <- "AC"
  cfs$tag_id_extra.1 <- "AC"

  out <- conform_genotypes(temp_csv(cfs), SYN_DWR, SYN_VCFS, verbose = FALSE,
                           keep_unmatched_markers = FALSE)

  expect_identical(colnames(out), colnames(syn_dwr()))

  txt <- paste(conform_report(out)$text, collapse = "\n")
  expect_match(txt, "did NOT map and will be dropped")

  expect_equal(conform_report(out)$coverage$extra_cols, character(0))
})

test_that("kept extra markers still get the same downstream QC", {
  # No VCF entry and no DWR column both mean no width to check against, so an
  # oddly-sized allele on a CFS-only marker is retained as-is: nothing claims
  # to have verified it, but nothing silently mangles it either.
  cfs <- syn_cfs()
  cfs$tag_id_extra <- "ACGT"
  cfs$tag_id_extra.1 <- "AC"

  out <- conform_genotypes(temp_csv(cfs), SYN_DWR, SYN_VCFS, verbose = FALSE)
  expect_true(all(out$tag_id_extra.1 == "ACGT"))
  expect_true(all(out$tag_id_extra.2 == "AC"))
})

test_that("a duplicate SAMPLE_ID's winner is chosen using extra markers too", {
  cfs <- syn_cfs()
  # The literal "NA" spelling, not NA_character_: temp_csv() round-trips
  # through write.csv(na = ""), which turns a real NA into an empty string on
  # read-back rather than back into NA, so a real NA here would not exercise
  # the "uncalled" path this test wants (see .is_called()'s doc comment for
  # why both spellings exist in real exports and must be treated alike).
  cfs$tag_id_extra <- "NA"
  cfs$tag_id_extra.1 <- "NA"

  # SYN_DUP_HIGHLOW's fuller row gets extra calls only a kept CFS-only marker
  # can provide; the emptier row stays empty. The winner should not change.
  rws <- which(cfs$indiv == "SYN_DUP_HIGHLOW")
  called <- rowSums(.is_called(cfs[rws, setdiff(colnames(cfs), SYN_META)]))
  fuller <- rws[which.max(called)]
  cfs$tag_id_extra[fuller] <- "AC"
  cfs$tag_id_extra.1[fuller] <- "AC"

  out <- conform_genotypes(temp_csv(cfs), SYN_DWR, SYN_VCFS, verbose = FALSE)
  kept <- out[out$SAMPLE_ID == "SYN_DUP_HIGHLOW", ]
  expect_equal(nrow(kept), 1L)
  expect_equal(kept$tag_id_extra.1, "AC")
})

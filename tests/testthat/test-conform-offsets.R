# The offset scan is the reason this conversion exists: a CFS column offset is
# only sometimes visible as a width change, so it has to be caught by allele
# validity instead. synthetic_cfs_offset.csv holds exactly the same genotype
# matrix as the clean fixture; only the header line is permuted.

test_that("the offset fixture differs from the clean one only in column order", {
  clean <- colnames(syn_cfs(SYN_CFS))
  broken <- colnames(syn_cfs(SYN_OFFSET))

  expect_setequal(clean, broken)
  expect_false(identical(clean, broken))
  # tag_id_101 was listed at the END of its block instead of the front, exactly
  # as tag_id_5720 was in the 20260313 export that motivated the scan.
  expect_lt(match("tag_id_101", clean), match("tag_id_101", broken))
  expect_lt(match("tag_id_110", broken), match("tag_id_101", broken))
})

test_that("the offset fixture is flagged, and the clean one is not", {
  broken <- syn_log(SYN_OFFSET)
  expect_match(broken, "samples carry a shifted block")
  expect_match(broken, "-1")
  expect_match(broken, "tag_id_10[3-9]")

  clean <- syn_log()
  expect_match(clean, "No shifted blocks detected")
  expect_false(grepl("samples carry a shifted block", clean))
})

test_that("the offset summary names the affected loci and their shift", {
  tab <- conform_report(syn_run(SYN_OFFSET))$offsets
  expect_s3_class(tab, "data.frame")
  # A contiguous run of 9 loci, all slipped one locus to the left, in most
  # samples. Sorted by extent, so the head is the worst-affected.
  expect_equal(nrow(tab), 9L)
  expect_true(all(tab$shift == -1L))
  expect_true(all(grepl("^tag_id_1", tab$locus)))
  expect_false(is.unsorted(rev(tab$samples)))

  expect_null(conform_report(syn_run())$offsets)
})

test_that("the mislabelled block is caught even where widths agree", {
  # tag_id_110 and tag_id_101 are both 2 wide, so the slip between them is
  # invisible to a length check and can only be seen by validating bases. If
  # this stops holding, the QC has quietly regressed to a width test.
  v <- syn_vcf()
  expect_equal(length(v[["tag_id_101"]]), length(v[["tag_id_110"]]))
  expect_match(syn_log(SYN_OFFSET), "tag_id_101")
})

test_that("an unrepaired offset is recorded as an anomaly", {
  expect_true(any(grepl("shifted genotype blocks and repair is OFF",
                        syn_anomalies(SYN_OFFSET))))
})

test_that("repair recovers genotypes the width filter would have discarded", {
  unrepaired <- syn_run(SYN_OFFSET)
  repaired <- syn_run(SYN_OFFSET, repair_offsets = TRUE)

  # Repair moves genotypes back to their own locus instead of letting the width
  # filter discard them, so it both retains far more calls and leaves far fewer
  # that are wrong for where they sit.
  expect_gt(n_called(repaired), n_called(unrepaired))

  bad_before <- count_invalid_alleles(unrepaired)
  bad_after <- count_invalid_alleles(repaired)
  expect_gt(bad_before, 0L)             # the fault is visible without repair
  expect_lt(bad_after, 0.5 * bad_before)

  # Some survive by design: the genotype shifted off the front of the block is
  # unrecoverable, and null_invalid_alleles is off, so nothing is silently
  # deleted. The caller can still see what is wrong.
  expect_gt(bad_after, 0L)
})

test_that("nulling invalid alleles is opt-in", {
  expect_gt(count_invalid_alleles(syn_run(SYN_OFFSET)), 0L)
  expect_equal(
    count_invalid_alleles(syn_run(SYN_OFFSET, repair_offsets = TRUE,
                                  null_invalid_alleles = TRUE)),
    0L
  )
})

test_that("the offset costs calls when it is left unrepaired", {
  # The point of the LOCI EMPTIED BY QC section: an unrepaired offset does not
  # look like corruption in the output, it looks like missing data, and the
  # report has to connect the two.
  txt <- syn_log(SYN_OFFSET)
  expect_match(txt, "LOCI EMPTIED BY QC")
  expect_match(txt, "repair_offsets = TRUE")
})

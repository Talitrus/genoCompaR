# What downstream consumers are entitled to assume about a conversion's output.

test_that("output columns match the DWR reference exactly, in order", {
  # The fixture stores DWR columns sorted and CFS columns in export order, so
  # this fails unless the conversion really reorders rather than just renames.
  expect_identical(colnames(syn_run()), colnames(syn_dwr()))
})

test_that("output satisfies the CKMRsim uniform-width invariant", {
  expect_equal(variable_width_loci(syn_run()), character(0))
})

test_that("every retained allele is a permitted allele of its locus", {
  expect_equal(count_invalid_alleles(syn_run()), 0L)
})

test_that("Ots_BMP-2-SNP1 carries data", {
  # The regression that motivated the locus: a rename rule replacing only the
  # first dot emits Ots_BMP-2.SNP1, which matches no DWR column, and the locus
  # comes out silently blank.
  out <- syn_run()
  expect_true(all(c("Ots_BMP-2-SNP1.1", "Ots_BMP-2-SNP1.2") %in% colnames(out)))
  expect_gt(sum(.is_called(out$`Ots_BMP-2-SNP1.1`)), 0)
})

test_that("DWR-only markers are present but empty", {
  out <- syn_run()
  renamed <- .transform_name(setdiff(colnames(syn_cfs()), SYN_META))
  dwr_only <- setdiff(setdiff(colnames(out), "SAMPLE_ID"), renamed)

  expect_gt(length(dwr_only), 0)
  # Greb1l / winter-run diagnostics are not on the CFS panel.
  expect_true(any(grepl("greb1l", dwr_only, ignore.case = TRUE)))
  for (col in dwr_only) {
    expect_true(all(!.is_called(out[[col]])),
                info = paste("expected empty:", col))
  }
})

test_that("a locus described by no VCF is declared, not silently trusted", {
  # tag_id_999 is in both CSVs and in neither VCF. Its width can only come from
  # the reference and its alleles cannot be checked at all; the report has to
  # say so rather than let it pass as verified.
  txt <- syn_log()
  expect_match(txt, "NO VCF entry")
  expect_match(txt, "alleles unchecked")
  expect_match(txt, "tag_id_999")

  coverage <- conform_report(syn_run())$coverage
  expect_equal(coverage$unverified, "tag_id_999")
  expect_length(coverage$validated, 19L)

  # It is still converted -- unverifiable is not the same as dropped.
  expect_gt(sum(.is_called(syn_run()$tag_id_999.1)), 0)
})

test_that("the report is retained whether or not it was printed", {
  quiet <- syn_run()
  expect_silent(conform_genotypes(SYN_CFS, SYN_DWR, SYN_VCFS, verbose = FALSE))
  expect_gt(length(conform_report(quiet)$text), 40)

  loud <- capture.output(
    out <- conform_genotypes(SYN_CFS, SYN_DWR, SYN_VCFS, verbose = TRUE)
  )
  expect_identical(loud, conform_report(out)$text)
})

test_that("conform_report() rejects anything else", {
  expect_error(conform_report(syn_dwr()), "must be a geno_conformed object")
})

test_that("the result prints as a conformed table, not a bare tibble", {
  out <- paste(capture.output(print(syn_run())), collapse = "\n")
  expect_match(out, "<geno_conformed> 37 samples x 45 columns \\| 1 anomaly")
  expect_match(out, "conform_report\\(\\)")
})

test_that("the result converts into the rest of the package", {
  geno <- as_geno_data(syn_run(), id_col = "SAMPLE_ID", group_col = NULL)
  expect_s3_class(geno, "geno_data")
  expect_equal(length(get_loci(geno)), 22L)
  expect_equal(nrow(geno), 37L)

  cmp <- compare_genotypes(geno, "SYN0001", "SYN0002", output = "tibble")
  expect_equal(cmp$n_loci_total[[1]], 22L)
})

test_that("call rates are reported per sample, named by ID", {
  rates <- conform_report(syn_run())$call_rates
  out <- syn_run()
  expect_named(rates, out$SAMPLE_ID)
  expect_equal(unname(rates[["SYN_NOCALL"]]), 0)
  expect_true(all(rates >= 0 & rates <= 1))
})

test_that("the written CSV round-trips to the same table", {
  out <- syn_run()
  tmp <- withr::local_tempfile(fileext = ".csv")
  write.csv(out, tmp, row.names = FALSE, na = "")
  back <- .read_chr(tmp)

  expect_identical(colnames(back), colnames(out))
  expect_equal(nrow(back), nrow(out))
  # write.csv(na = "") then read.csv turns empties back into NA; compare calls.
  expect_equal(sum(.is_called(unlist(back[, -1], use.names = FALSE))),
               sum(.is_called(unlist(out[, -1], use.names = FALSE))))
})

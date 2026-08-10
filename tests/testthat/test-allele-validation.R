test_that(".is_called() treats NA, empty and literal 'NA' as uncalled", {
  x <- c("AG", NA, "", "NA", "T")
  expect_equal(.is_called(x), c(TRUE, FALSE, FALSE, FALSE, TRUE))
})

test_that(".valid_against() accepts only permitted bases at each site", {
  allow <- list(c("A", "G"), c("C", "T"))
  expect_equal(.valid_against(c("AC", "GT", "AT", "GC"), allow), rep(TRUE, 4))
  # Wrong base at a site.
  expect_equal(.valid_against(c("TC", "AA"), allow), c(FALSE, FALSE))
})

test_that(".valid_against() rejects the wrong width", {
  allow <- list(c("A", "G"), c("C", "T"))
  expect_equal(.valid_against(c("A", "ACT", "AC"), allow), c(FALSE, FALSE, TRUE))
})

test_that(".valid_against() reports uncalled values as invalid, not NA", {
  allow <- list(c("A", "G"))
  res <- .valid_against(c(NA, "", "NA", "A"), allow)
  expect_equal(res, c(FALSE, FALSE, FALSE, TRUE))
  expect_false(anyNA(res))
})

test_that(".valid_against() handles an all-uncalled vector", {
  allow <- list(c("A", "G"))
  expect_equal(.valid_against(c(NA, NA), allow), c(FALSE, FALSE))
})

test_that(".valid_allele() falls back on alleles seen in the reference", {
  permitted <- list(allow = list(c("A", "G")), observed = c("T"))
  expect_equal(.valid_allele(c("A", "T", "C"), permitted), c(TRUE, TRUE, FALSE))
  # A locus the VCF does not describe permits nothing.
  expect_equal(.valid_allele(c("A", "T"), NULL), c(FALSE, FALSE))
})

test_that("reference allele widths match the VCF site counts", {
  # The assumption the whole conversion rests on: the VCFs and the reference
  # agree, so when CFS disagrees with them it is CFS that is wrong.
  v <- syn_vcf()
  dwr <- syn_dwr()
  mismatched <- character(0)
  for (locus in names(v)) {
    col <- paste0(locus, ".1")
    if (!col %in% colnames(dwr)) next
    vals <- dwr[[col]]
    vals <- vals[.is_called(vals)]
    if (length(vals) == 0) next
    widths <- unique(nchar(vals))
    if (length(widths) != 1 || widths != length(v[[locus]])) {
      mismatched <- c(mismatched, locus)
    }
  }
  expect_equal(mismatched, character(0),
               info = paste("reference disagrees with VCF at:",
                            paste(head(mismatched, 10), collapse = ", ")))
})

test_that("the VCF ALT lists nearly, but not fully, cover the reference alleles", {
  # The cohort-called-VCF problem, reproduced: one real reference allele
  # (tag_id_106 "TCAGA") is absent from every ALT list. The gap is small and
  # mixed in origin - some genuinely rare variants, some miscalls, and
  # frequency alone does not separate them. .build_permitted() therefore widens
  # rather than hard-filters, and null_invalid_alleles defaults to FALSE, so a
  # rare real allele is never deleted as corruption.
  v <- syn_vcf()
  dwr <- syn_dwr()
  bad <- 0L
  total <- 0L
  for (locus in names(v)) {
    for (col in paste0(locus, c(".1", ".2"))) {
      if (!col %in% colnames(dwr)) next
      vals <- dwr[[col]]
      vals <- vals[.is_called(vals)]
      if (length(vals) == 0) next
      total <- total + length(vals)
      bad <- bad + sum(!.valid_against(vals, v[[locus]]))
    }
  }
  expect_gt(total, 0)
  expect_gt(bad, 0L)              # documents the gap rather than hiding it
  expect_lt(bad / total, 0.01)    # but it is small enough to rely on
})

test_that("widened permitted sets accept every reference allele", {
  expect_equal(count_invalid_alleles(syn_dwr()), 0L)
})

test_that("widening never admits an allele of the wrong width", {
  permitted <- syn_permitted()
  for (locus in names(permitted)) {
    obs <- permitted[[locus]]$observed
    if (length(obs) == 0) next
    expect_true(all(nchar(obs) == length(permitted[[locus]]$allow)),
                info = paste("wrong-width allele admitted at", locus))
  }
})

test_that("the reference satisfies the CKMRsim uniform-width invariant", {
  expect_equal(variable_width_loci(syn_dwr()), character(0))
})

test_that("the clean CFS export has no unrecognised alleles", {
  # A canary on incoming data: an allele that is neither VCF-permitted nor seen
  # anywhere in the reference is either novel or a corruption, and either way
  # wants a human to look before the file is trusted.
  permitted <- syn_permitted()
  cfs <- syn_cfs()
  bad <- 0L
  total <- 0L
  for (raw in setdiff(colnames(cfs), SYN_META)) {
    locus <- sub("\\.[12]$", "", .transform_name(raw))
    if (is.null(permitted[[locus]])) next
    vals <- cfs[[raw]]
    vals <- vals[.is_called(vals)]
    if (length(vals) == 0) next
    total <- total + length(vals)
    bad <- bad + sum(!.valid_allele(vals, permitted[[locus]]))
  }
  expect_gt(total, 0)
  expect_equal(bad, 0L,
               info = sprintf("%d of %d CFS alleles are unrecognised", bad, total))
})

test_that("the fixture carries both spellings of an uncalled genotype", {
  # Exports have used the literal string NA and an empty cell interchangeably;
  # .is_called() has to treat them alike, so the fixture contains both.
  raw <- readLines(SYN_CFS)
  expect_true(any(grepl(",NA,", raw)))
  expect_true(any(grepl(",,", raw)))
})

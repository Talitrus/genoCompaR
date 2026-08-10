test_that("the reference VCFs load into per-locus allele definitions", {
  v <- syn_vcf()
  expect_type(v, "list")
  expect_equal(length(v), 19L)
  # Each locus is a list of per-site character vectors of permitted bases.
  expect_type(v[[1]], "list")
  expect_type(v[[1]][[1]], "character")

  # Every width the real panel spans, in miniature: a shift between two loci of
  # equal width cannot be caught by length, which is why widths alone are not
  # the QC.
  expect_setequal(unique(vapply(v, length, integer(1))), c(1L, 2L, 3L, 4L, 5L))
})

test_that("loci repeated across the two VCFs are not double counted", {
  # The panel VCF repeats loci already present in the full VCF at identical
  # (CHROM, POS). Counting both copies doubles the apparent haplotype width and
  # makes every real allele look wrong, so dedup is load-bearing.
  v <- syn_vcf()
  shared <- c("NC_037099.1:62937268-62937373", "NC_037130.1:1062935-1063235")
  expect_true(all(shared %in% names(v)))

  only_full <- .load_vcf_alleles(SYN_VCFS[1])
  for (locus in shared) {
    expect_equal(length(v[[locus]]), length(only_full[[locus]]),
                 info = paste("width changed after merging VCFs at", locus))
    expect_equal(length(v[[locus]]), 2L)   # 4 would mean the dedup was lost
  }
})

test_that("indel records never become haplotype sites", {
  v <- syn_vcf()

  # This locus carries two SNP sites and two indel records. Counting the indels
  # would make it 4 wide and let multi-base values in.
  expect_equal(length(v[["NC_037124.1:12281207-12281551"]]), 2L)

  # This one is nothing but indels, so after filtering it must be gone rather
  # than present and zero wide.
  expect_false("NC_037124.1:12270118-12270418" %in% names(v))

  for (locus in names(v)) {
    for (site in v[[locus]]) {
      expect_true(all(nchar(site) == 1L),
                  info = paste("multi-base base at", locus))
    }
  }
})

test_that("a multi-base ALT at a SNP site is not a permitted base", {
  # tag_id_104's second site is REF A, ALT "G,AT". The site stays; "AT" does not
  # become a value the haplotype may take there.
  expect_setequal(syn_vcf()[["tag_id_104"]][[2]], c("A", "G"))
})

test_that("every site offers at least one permitted base", {
  for (locus in syn_vcf()) {
    for (site in locus) expect_gte(length(site), 1)
  }
})

test_that("a missing VCF degrades coverage rather than raising", {
  # Losing one reference still leaves the others validating what they cover,
  # but it silently shrinks QC coverage, so it has to be surfaced.
  expect_warning(
    v <- .load_vcf_alleles(c(SYN_VCFS[1], syn_path("no_such.vcf"))),
    "VCF not found"
  )
  expect_gt(length(v), 0)
})

test_that("no usable VCF at all is fatal", {
  expect_error(
    suppressWarnings(.load_vcf_alleles(syn_path("no_such.vcf"))),
    "No VCF records loaded"
  )
})

test_that("the shipped reference VCFs describe the full Chinook panel", {
  # These are what conform_genotypes() is run against in production.
  vcfs <- system.file(
    "extdata", "reference_VCFs",
    c("full.vcf", "Chinook_FullPanel_VGLL3Six6LFARWRAP.vcf"),
    package = "genoCompaR"
  )
  expect_length(vcfs, 2L)

  v <- .load_vcf_alleles(vcfs)
  expect_equal(length(v), 203L)
  # Spot checks against the VCF records themselves.
  expect_equal(length(v[["tag_id_684"]]), 2L)
  expect_equal(length(v[["tag_id_871"]]), 9L)
  expect_equal(length(v[["NC_037130.1:1062935-1063235"]]), 2L)
})

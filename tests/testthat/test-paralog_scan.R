# Heterozygote excess catches a class of miscall nothing else here can see: two
# paralogous loci collapsed into one call. Every genotype is individually
# well-formed - right width, permitted bases - but nearly every fish comes out
# heterozygous, which no allele frequency distribution permits.

one_locus <- function(a1, a2) {
  data.frame(L.1 = a1, L.2 = a2, stringsAsFactors = FALSE)
}

test_that("expected heterozygosity is bounded by 0.5 for a biallelic locus", {
  # Sanity check on the statistic itself before trusting it on real data.
  s <- paralog_scan(one_locus(c("A", "A", "G", "G"), c("A", "G", "G", "A")),
                    min_n = 1)
  expect_equal(s$exp_het, 0.5)
  expect_equal(s$obs_het, 0.5)
  expect_equal(s$alleles, 2L)
})

test_that("a fully heterozygous biallelic locus is impossible under HWE", {
  s <- paralog_scan(one_locus(rep("A", 100), rep("G", 100)))
  expect_equal(s$obs_het, 1)
  expect_equal(s$exp_het, 0.5)
  expect_equal(s$excess, 0.5)
  expect_true(s$flagged)
})

test_that("a multi-allelic microhap may legitimately exceed 0.5 heterozygosity", {
  # Guards against a naive "het > 0.5 is broken" rule, which would condemn
  # perfectly good microhaplotype loci. Only the gap from the locus's own
  # expectation is evidence of a problem.
  alleles <- c("AAAA", "CCCC", "GGGG", "TTTT", "ACGT", "TGCA")
  set.seed(1)
  s <- paralog_scan(one_locus(sample(alleles, 400, replace = TRUE),
                              sample(alleles, 400, replace = TRUE)))
  expect_gt(s$obs_het, 0.5)
  expect_false(s$flagged)
})

test_that("loci with too few genotypes are not judged", {
  # Small n makes the statistic meaningless; the check must skip rather than
  # emit noise.
  s <- paralog_scan(one_locus(rep("A", 5), rep("G", 5)))
  expect_equal(nrow(s), 0L)
  expect_equal(attr(s, "n_loci"), 1L)

  expect_equal(nrow(paralog_scan(one_locus(rep("A", 5), rep("G", 5)), min_n = 5)), 1L)
})

test_that("uncalled genotypes are excluded, in every spelling", {
  s <- paralog_scan(one_locus(c("A", "A", NA, "G", "A"), c("G", "G", "G", "", "NA")),
                    min_n = 1)
  expect_equal(s$n, 2L)
  expect_equal(s$obs_het, 1)
})

test_that("results are sorted by excess, worst first", {
  df <- data.frame(
    calm.1 = rep(c("A", "A", "G", "G"), 25),
    calm.2 = rep(c("A", "G", "G", "A"), 25),
    collapsed.1 = rep("A", 100),
    collapsed.2 = rep("G", 100),
    stringsAsFactors = FALSE
  )
  s <- paralog_scan(df)
  expect_equal(s$locus, c("collapsed", "calm"))
  expect_true(!is.unsorted(rev(s$excess)))
  expect_equal(sum(s$flagged), 1L)
})

test_that("the threshold controls what is flagged, not what is reported", {
  df <- data.frame(a.1 = rep("A", 100), a.2 = rep("G", 100),
                   stringsAsFactors = FALSE)
  expect_true(paralog_scan(df, threshold = 0.25)$flagged)
  expect_false(paralog_scan(df, threshold = 0.9)$flagged)
  expect_equal(nrow(paralog_scan(df, threshold = 0.9)), 1L)
})

test_that("a geno_data, a conformed table and a bare data frame all work", {
  geno <- synthetic_geno()
  s <- paralog_scan(geno, min_n = 1)
  expect_s3_class(s, "geno_paralog_scan")
  expect_setequal(s$locus, get_loci(geno))

  # A bare data frame's loci are derived from its column names, so the ID and
  # group columns must not be mistaken for one.
  s2 <- paralog_scan(as.data.frame(geno), min_n = 1)
  expect_setequal(s2$locus, get_loci(geno))

  # A conformed table is neither, and its locus names contain their own dots.
  s3 <- paralog_scan(syn_run(), min_n = 1)
  expect_true("NC_037099.1:62937268-62937373" %in% s3$locus)
  expect_false(any(grepl("SAMPLE_ID", s3$locus)))
})

test_that("an explicit locus subset is honoured", {
  geno <- synthetic_geno()
  s <- paralog_scan(geno, min_n = 1, loci = c("locus1", "locus2"))
  expect_setequal(s$locus, c("locus1", "locus2"))
  expect_equal(attr(s, "n_loci"), 2L)
})

test_that("input without allele-column pairs is an error, not an empty result", {
  expect_error(paralog_scan(data.frame(a = 1, b = 2)), "No locus allele-column")
  expect_error(paralog_scan("not a data frame"), "must be a data frame")
})

test_that("the print method reports what was scanned and what was skipped", {
  df <- data.frame(
    big.1 = rep("A", 100), big.2 = rep("G", 100),
    small.1 = c(rep("A", 5), rep(NA, 95)), small.2 = c(rep("G", 5), rep(NA, 95)),
    stringsAsFactors = FALSE
  )
  out <- paste(capture.output(print(paralog_scan(df))), collapse = "\n")
  expect_match(out, "1 of 2 loci scanned, 1 flagged")
  expect_match(out, "1 locus/loci skipped")
})

test_that("the broodstock reference carries two collapsed-paralog loci", {
  # Both are called at ~93-100% heterozygosity from two alleles, which cannot
  # happen. Recorded here so a future reference build that fixes - or spreads -
  # the problem is noticed.
  ref <- .read_chr(system.file(
    "extdata", "Parents_FRH2024_broodstock_reference_genotypes.csv",
    package = "genoCompaR"
  ))
  s <- paralog_scan(ref)

  expect_setequal(s$locus[s$flagged],
                  c("NC_037104.1:56552952-56553042",
                    "NC_037130.1:1062935-1063235"))
  flagged <- s[s$flagged, ]
  expect_true(all(flagged$alleles == 2L))
  expect_true(all(flagged$obs_het > 0.9))
  expect_true(all(flagged$n > 1000L))
})

test_that("conform_genotypes() reports heterozygote excess as an anomaly", {
  # The conversion runs the same scan over its own output. The clean fixture
  # carries no collapsed paralog, so one is injected here: every sample made
  # heterozygous for the same two valid alleles at tag_id_101. Both alleles are
  # taken from the reference, so the width and validity checks pass them and
  # only the heterozygosity check can object.
  dwr <- syn_dwr()
  alleles <- unique(c(dwr$tag_id_101.1, dwr$tag_id_101.2))
  alleles <- alleles[.is_called(alleles)]
  expect_gte(length(alleles), 2L)

  cfs <- syn_cfs()
  cfs$tag_id_101 <- alleles[1]      # CFS spells the pair "locus" / "locus.1"
  cfs$tag_id_101.1 <- alleles[2]

  out <- conform_genotypes(temp_csv(cfs), SYN_DWR, SYN_VCFS,
                           min_het_samples = 20L, verbose = FALSE)
  txt <- paste(conform_report(out)$text, collapse = "\n")

  expect_match(txt, "HETEROZYGOTE EXCESS")
  expect_match(txt, "collapsed paralogs")
  expect_match(txt, "tag_id_101")
  expect_true(any(grepl("heterozygote excess", conform_report(out)$anomalies)))

  # Nothing else objected: the calls themselves are well-formed.
  expect_equal(count_invalid_alleles(out), 0L)
  expect_equal(variable_width_loci(out), character(0))

  hit <- conform_report(out)$het_excess
  hit <- hit[hit$locus == "tag_id_101", ]
  expect_equal(hit$obs_het, 1)
  expect_true(hit$flagged)
})

test_that("the report always carries the full scan, flagged or not", {
  het <- conform_report(syn_run())$het_excess
  expect_s3_class(het, "geno_paralog_scan")
  expect_equal(attr(het, "threshold"), 0.25)
  expect_equal(as.integer(attr(het, "min_n")), 50L)
})

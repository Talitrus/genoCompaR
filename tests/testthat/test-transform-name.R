test_that("allele suffixes shift from CFS to DWR convention", {
  # CFS writes "locus" / "locus.1"; DWR writes "locus.1" / "locus.2".
  expect_equal(.transform_name("tag_id_1030"), "tag_id_1030.1")
  expect_equal(.transform_name("tag_id_1030.1"), "tag_id_1030.2")
})

test_that("NC_ genomic coordinates become colon/dash form", {
  expect_equal(.transform_name("NC_037099.1.62937268.62937373"),
               "NC_037099.1:62937268-62937373.1")
  expect_equal(.transform_name("NC_037099.1.62937268.62937373.1"),
               "NC_037099.1:62937268-62937373.2")
  # The assembly version (".1" after the accession) must survive intact.
  expect_true(grepl("^NC_037104\\.1:",
                    .transform_name("NC_037104.1.55923357.55923657")))
})

test_that("Ots_ and OkiOts_ markers dash throughout", {
  expect_equal(.transform_name("Ots_111312.435"), "Ots_111312-435.1")
  expect_equal(.transform_name("Ots_111312.435.1"), "Ots_111312-435.2")
  expect_equal(.transform_name("OkiOts_120255.113"), "OkiOts_120255-113.1")
})

test_that("Ots_BMP.2.SNP1 dashes every dot", {
  # Regression: an earlier version of this rule replaced only the first dot,
  # producing "Ots_BMP-2.SNP1", which matches nothing in the DWR panel. The
  # locus was then silently blank in its output.
  expect_equal(.transform_name("Ots_BMP.2.SNP1"), "Ots_BMP-2-SNP1.1")
  expect_equal(.transform_name("Ots_BMP.2.SNP1.1"), "Ots_BMP-2-SNP1.2")
})

test_that("u07 markers keep their trailing dot", {
  # DWR spells these "Ots_u07-07.161", not "Ots_u07-07-161", so u07 is the one
  # family where only the first dot becomes a dash.
  expect_equal(.transform_name("Ots_u07.07.161"), "Ots_u07-07.161.1")
  expect_equal(.transform_name("Ots_u07.07.161.1"), "Ots_u07-07.161.2")
  expect_equal(.transform_name("Ots_u07.49.290"), "Ots_u07-49.290.1")
})

test_that("names outside the known families pass through unchanged", {
  expect_equal(.transform_name("some_new_marker"), "some_new_marker.1")
})

test_that("the rule is vectorized and unnamed", {
  out <- .transform_name(c("tag_id_1030", "Ots_111312.435.1"))
  expect_equal(out, c("tag_id_1030.1", "Ots_111312-435.2"))
  expect_null(names(out))
})

test_that("the fixture covers every marker-naming family", {
  cfs_cols <- setdiff(colnames(syn_cfs()), SYN_META)

  # One representative of each rule .transform_name() implements.
  expect_true("NC_037099.1.62937268.62937373" %in% cfs_cols)  # coords
  expect_true("OkiOts_120255.113" %in% cfs_cols)              # OkiOts_
  expect_true("Ots_111312.435" %in% cfs_cols)                 # Ots_
  expect_true("Ots_BMP.2.SNP1" %in% cfs_cols)                 # every dot dashed
  expect_true("Ots_u07.07.161" %in% cfs_cols)                 # u07 exception
  expect_true("tag_id_101" %in% cfs_cols)                     # suffix only
})

test_that("every CFS marker in the fixture maps onto a DWR column", {
  renamed <- .transform_name(setdiff(colnames(syn_cfs()), SYN_META))

  unmapped <- setdiff(renamed, colnames(syn_dwr()))
  expect_equal(unmapped, character(0),
               info = paste("unmapped:", paste(head(unmapped, 10), collapse = ", ")))
  expect_false(any(duplicated(renamed)))
})

SYN_NMFS_PBT <- system.file(
  "extdata", "synthetic", "synthetic_nmfs_pbt.xlsx",
  package = "genoCompaR"
)

EXPECTED_PO_COLS <- c(
  "D2_indiv",
  "D1_indiv",
  "logl_ratio",
  "num_loc",
  "total_incompat",
  "total_compat",
  "fraction_incompat"
)

test_that("nmfs_pbt2ckmrout validates arguments", {
  expect_error(
    nmfs_pbt2ckmrout("nonexistent_file_12345.xlsx"),
    "File not found"
  )
  expect_error(
    nmfs_pbt2ckmrout(123),
    "must be a path to an Excel file or a data.frame"
  )
  expect_error(
    nmfs_pbt2ckmrout(data.frame(ma = "M1", kid1 = "K1")),
    "must contain 'ma' and 'pa'"
  )
  expect_error(
    nmfs_pbt2ckmrout(data.frame(ma = "M1", pa = "P1")),
    "at least one 'kid' column"
  )
  expect_error(
    nmfs_pbt2ckmrout(
      data.frame(ma = "M1", pa = "P1", kid1 = "K1"),
      strip_kid_suffix = NA
    ),
    "strip_kid_suffix"
  )
  expect_error(
    nmfs_pbt2ckmrout(
      data.frame(ma = "M1", pa = "P1", kid1 = "K1"),
      out_file = ""
    ),
    "out_file"
  )
})

test_that("nmfs_pbt2ckmrout expands single-parent, two-parent, and multi-kid families", {
  df <- data.frame(
    index = c(1, 2, 3, NA),
    ma = c("C250002MER", "NA", "T166266_7949_scarfF_2018", NA),
    pa = c("NA", "T157393_5246_scarfM_2017", "T246584_SJRRP2025ADULT132_returnM", NA),
    fs_fam_size = c(2, 1, 2, 5),
    kid1 = c("C260698CVP_LNNLLNLLLLLL", "C251376CVP_ENNEENEEEEEE", "C259512CVP_ENNEENEEEEEE", NA),
    kid2 = c("C260734CVP_LNNLLNLLLLLL", NA, "C259513CVP_EHHEEHEEEEEH", NA),
    kid3 = c(NA, NA, NA, NA),
    stringsAsFactors = FALSE
  )

  res <- nmfs_pbt2ckmrout(df)

  expect_s3_class(res, "tbl_df")
  expect_equal(names(res), EXPECTED_PO_COLS)
  # Family 1: 2 kids x 1 ma = 2 pairs
  # Family 2: 1 kid x 1 pa = 1 pair
  # Family 3: 2 kids x (ma + pa) = 4 pairs
  # Row 4 (total row): ignored
  expect_equal(nrow(res), 7L)

  expect_equal(
    res$D2_indiv,
    c(
      "C260698CVP", "C260734CVP",
      "C251376CVP",
      "C259512CVP", "C259512CVP",
      "C259513CVP", "C259513CVP"
    )
  )
  expect_equal(
    res$D1_indiv,
    c(
      "C250002MER", "C250002MER",
      "T157393_5246_scarfM_2017",
      "T166266_7949_scarfF_2018", "T246584_SJRRP2025ADULT132_returnM",
      "T166266_7949_scarfF_2018", "T246584_SJRRP2025ADULT132_returnM"
    )
  )

  # Missing CKMRsim columns are NA with appropriate types
  expect_type(res$logl_ratio, "double")
  expect_type(res$num_loc, "integer")
  expect_type(res$total_incompat, "integer")
  expect_type(res$total_compat, "integer")
  expect_type(res$fraction_incompat, "double")
  expect_true(all(is.na(res$logl_ratio)))
  expect_true(all(is.na(res$num_loc)))
  expect_true(all(is.na(res$total_incompat)))
  expect_true(all(is.na(res$total_compat)))
  expect_true(all(is.na(res$fraction_incompat)))
})

test_that("strip_kid_suffix = FALSE preserves full kid IDs", {
  df <- data.frame(
    ma = "C250002MER",
    pa = "NA",
    kid1 = "C260698CVP_LNNLLNLLLLLL",
    stringsAsFactors = FALSE
  )

  res_keep <- nmfs_pbt2ckmrout(df, strip_kid_suffix = FALSE)
  expect_equal(res_keep$D2_indiv, "C260698CVP_LNNLLNLLLLLL")
  expect_equal(res_keep$D1_indiv, "C250002MER")
})

test_that("empty or summary-only input returns 0-row tibble with expected schema", {
  df <- data.frame(
    index = NA,
    ma = NA,
    pa = "NA",
    fs_fam_size = 90,
    kid1 = NA,
    stringsAsFactors = FALSE
  )

  res <- nmfs_pbt2ckmrout(df)
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 0L)
  expect_equal(names(res), EXPECTED_PO_COLS)
  expect_type(res$D2_indiv, "character")
  expect_type(res$D1_indiv, "character")
  expect_type(res$logl_ratio, "double")
  expect_type(res$num_loc, "integer")
})

test_that("out_file writes a tab-delimited file matching Coleman_BY25_PO_results.tsv format", {
  df <- data.frame(
    ma = c("C250002MER", "NA"),
    pa = c("NA", "25-11669"),
    kid1 = c("C260698CVP_LNNLLNLLLLLL", "C260291SWP_ENNEENEEEEEE"),
    stringsAsFactors = FALSE
  )

  tmp_tsv <- withr::local_tempfile(fileext = ".tsv")
  res <- nmfs_pbt2ckmrout(df, out_file = tmp_tsv)

  expect_true(file.exists(tmp_tsv))
  lines <- readLines(tmp_tsv)
  expect_equal(
    lines[1],
    paste(EXPECTED_PO_COLS, collapse = "\t")
  )
  expect_equal(
    lines[2],
    "C260698CVP\tC250002MER\tNA\tNA\tNA\tNA\tNA"
  )
  expect_equal(
    lines[3],
    "C260291SWP\t25-11669\tNA\tNA\tNA\tNA\tNA"
  )

  parsed <- utils::read.delim(tmp_tsv, stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(names(parsed), EXPECTED_PO_COLS)
  expect_equal(parsed$D2_indiv, res$D2_indiv)
  expect_equal(parsed$D1_indiv, res$D1_indiv)
})

test_that("nmfs_pbt2ckmrout reads Excel workbooks (.xlsx) and skips summary rows", {
  xlsx_path <- if (nzchar(SYN_NMFS_PBT) && file.exists(SYN_NMFS_PBT)) {
    SYN_NMFS_PBT
  } else {
    testthat::test_path("..", "..", "inst", "extdata", "synthetic", "synthetic_nmfs_pbt.xlsx")
  }
  skip_if_not(file.exists(xlsx_path))

  res <- nmfs_pbt2ckmrout(xlsx_path)

  expect_s3_class(res, "tbl_df")
  expect_equal(names(res), EXPECTED_PO_COLS)
  # Row 2: 2 kids x 1 ma = 2
  # Row 3: 1 kid x 1 pa = 1
  # Row 4: 2 kids x 2 parents = 4
  # Row 5: summary row (fs_fam_size = 5) skipped
  expect_equal(nrow(res), 7L)
  expect_equal(
    res$D2_indiv,
    c("KID_01", "KID_02", "KID_03", "KID_04", "KID_04", "KID_05", "KID_05")
  )
  expect_equal(
    res$D1_indiv,
    c(
      "MA_001", "MA_001",
      "PA_002_scarfM_2018",
      "MA_003_scarfF_2019", "PA_003_returnM",
      "MA_003_scarfF_2019", "PA_003_returnM"
    )
  )

  # Reading a non-family sheet errors cleanly
  expect_error(
    nmfs_pbt2ckmrout(xlsx_path, sheet = "unassigned"),
    "must contain 'ma' and 'pa'"
  )
})

test_that("nmfs_pbt2ckmrout converts example_data/Salvage2026_SanJoaquin_ckmr.xlsx when present", {
  ex_xlsx <- testthat::test_path("..", "..", "example_data", "Salvage2026_SanJoaquin_ckmr.xlsx")
  ex_tsv <- testthat::test_path("..", "..", "example_data", "Coleman_BY25_PO_results.tsv")
  skip_if_not(file.exists(ex_xlsx))

  res <- nmfs_pbt2ckmrout(ex_xlsx)

  expect_equal(nrow(res), 105L)
  expect_equal(res$D2_indiv[1:2], c("C260698CVP", "C260734CVP"))
  expect_equal(res$D1_indiv[1:2], c("C250002MER", "C250002MER"))
  expect_false(any(grepl("_", res$D2_indiv)))
  expect_true(all(is.na(res$logl_ratio)))

  if (file.exists(ex_tsv)) {
    coleman <- utils::read.delim(ex_tsv, stringsAsFactors = FALSE, check.names = FALSE)
    expect_equal(names(res), names(coleman))
  }
})


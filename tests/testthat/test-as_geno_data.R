test_that("as_geno_data classes and validates a plain data frame", {
  df <- synthetic_geno_df()
  g <- as_geno_data(df)
  expect_s3_class(g, "geno_data")
  expect_equal(get_loci(g), paste0("locus", 1:5))
  expect_equal(nrow(g), 6)
})

test_that("as_geno_data can restrict to a locus subset", {
  g <- as_geno_data(synthetic_geno_df(), loci = c("locus1", "locus3"))
  expect_equal(get_loci(g), c("locus1", "locus3"))
  expect_false("locus2.1" %in% names(g))
})

test_that("as_geno_data rejects duplicate IDs and missing columns", {
  df <- synthetic_geno_df()
  df$indiv.ID[2] <- df$indiv.ID[1]
  expect_error(as_geno_data(df), "Duplicate sample IDs")
  expect_error(as_geno_data(synthetic_geno_df(), id_col = "nope"), "not found")
  expect_error(as_geno_data(synthetic_geno_df(), loci = "nosuch"), "missing allele column")
  expect_error(as_geno_data("not a data frame"), "data frame")
})

test_that("merge_genotypes combines two sets on shared loci", {
  a <- as_geno_data(synthetic_geno_df()[1:3, ], loci = paste0("locus", 1:5))
  b <- as_geno_data(synthetic_geno_df()[4:6, ], loci = paste0("locus", 1:3))

  m <- merge_genotypes(a, b)
  expect_s3_class(m, "geno_data")
  expect_equal(nrow(m), 6)
  expect_equal(get_loci(m), paste0("locus", 1:3)) # intersection only
  expect_setequal(m$indiv.ID, synthetic_geno_df()$indiv.ID)
})

test_that("merge_genotypes reports how many loci are unique to each side", {
  a <- as_geno_data(synthetic_geno_df()[1:3, ], loci = paste0("locus", 1:5))
  b <- as_geno_data(synthetic_geno_df()[4:6, ], loci = paste0("locus", 1:3))

  expect_message(merge_genotypes(a, b),
                 "2 locus/loci only in `x` and 0 only in `y`.*dropped")
  expect_message(merge_genotypes(a, b, keep_unmatched_loci = TRUE),
                 "2 locus/loci only in `x` and 0 only in `y`.*kept, NA-filled")

  # A clean intersection (no unique loci on either side) has nothing to
  # report.
  c <- as_geno_data(synthetic_geno_df()[4:6, ], loci = paste0("locus", 1:5))
  expect_silent(merge_genotypes(a, c))
})

test_that("keep_unmatched_loci keeps the union, NA-filled where absent", {
  a <- as_geno_data(synthetic_geno_df()[1:3, ], loci = paste0("locus", 1:5))
  b <- as_geno_data(synthetic_geno_df()[4:6, ], loci = paste0("locus", 1:3))

  m <- suppressMessages(merge_genotypes(a, b, keep_unmatched_loci = TRUE))
  expect_equal(get_loci(m), paste0("locus", 1:5))
  expect_equal(nrow(m), 6)

  b_rows <- m[m$indiv.ID %in% c("unrelated1", "dup_a", "dup_b"), ]
  expect_true(all(is.na(b_rows$locus4.1) & is.na(b_rows$locus4.2)))
  expect_true(all(is.na(b_rows$locus5.1) & is.na(b_rows$locus5.2)))

  a_rows <- m[m$indiv.ID %in% c("mom", "dad", "kid"), ]
  expect_false(any(is.na(a_rows$locus3.1)))  # a shared locus is unaffected

  # A locus unique to one side is simply excluded from any comparison that
  # crosses the merge - not an error, not a false incompatibility.
  cmp <- compare_genotypes(m, "kid", "unrelated1", output = "tibble")
  expect_equal(cmp$n_loci_total[[1]], 5L)
  expect_equal(cmp$n_loci_compared[[1]], 3L)
})

test_that("keep_unmatched_loci widens what an explicit loci argument accepts", {
  a <- as_geno_data(synthetic_geno_df()[1:3, ], loci = paste0("locus", 1:5))
  b <- as_geno_data(synthetic_geno_df()[4:6, ], loci = paste0("locus", 1:3))

  expect_error(merge_genotypes(a, b, loci = "locus4"),
               "not present in both inputs")

  m <- merge_genotypes(a, b, loci = "locus4", keep_unmatched_loci = TRUE)
  expect_equal(get_loci(m), "locus4")
  expect_error(
    merge_genotypes(a, b, loci = "nosuchlocus", keep_unmatched_loci = TRUE),
    "not present in either input"
  )
})

test_that("merge_genotypes rejects overlapping IDs and disjoint panels", {
  a <- as_geno_data(synthetic_geno_df()[1:3, ])
  expect_error(merge_genotypes(a, a), "Duplicate sample IDs")

  df <- synthetic_geno_df()
  names(df) <- sub("^locus", "other", names(df))
  b <- as_geno_data(df, loci = paste0("other", 1:5))
  expect_error(merge_genotypes(a, b), "share no loci")
  expect_error(merge_genotypes(a, unclass(a)), "geno_data")
})

test_that("merged genotypes reproduce the known real trio result", {
  # the API-level replacement for hand-assembling a geno_data: read each
  # file, subset, merge - and get the same answer as the manual fixture
  geno_manual <- .known_trio_geno()

  parents <- read_genotypes(
    .parents_file_path(), id_col = "SAMPLE_ID", group_col = NULL
  )
  raw <- utils::read.csv(
    .unknowns_file_path(), colClasses = "character",
    na.strings = c("NA", ""), check.names = FALSE, stringsAsFactors = FALSE
  )
  offspring <- as_geno_data(
    raw[raw$SAMPLE_ID == "C260383CVP", ],
    id_col = "SAMPLE_ID", group_col = NULL
  )
  merged <- merge_genotypes(
    offspring,
    parents[parents$SAMPLE_ID %in% c("C240716FRH", "C240722FRH"), ]
  )

  expect_equal(nrow(merged), 3)
  expect_equal(length(get_loci(merged)), length(get_loci(geno_manual)))

  a <- check_trio(merged, "C260383CVP", "C240716FRH", "C240722FRH", output = "tibble")
  b <- check_trio(geno_manual, "C260383CVP", "C240716FRH", "C240722FRH", output = "tibble")
  expect_equal(a$n_compatible, b$n_compatible)
  expect_equal(a$n_incompatible, b$n_incompatible)
  expect_equal(a$n_loci_compared, b$n_loci_compared)
})

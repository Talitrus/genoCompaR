# CFS re-runs a failed sample under its original ID, so an ID can appear more
# than once. Agreed rule: keep the row with the most calls; drop the ID
# entirely when every one of its rows is empty.

test_that("the fixture contains all four duplicate-ID shapes", {
  ids <- syn_cfs()[[SYN_ID_COL]]
  expect_equal(sum(ids == "SYN_DUP_HIGHLOW"), 2L)   # failed run + good run
  expect_equal(sum(ids == "SYN_DUP_EMPTY"), 2L)     # both runs failed
  expect_equal(sum(ids == "SYN_DUP_DISCORD"), 2L)   # two good runs that disagree
  expect_equal(sum(ids == "SYN_DUP_TRIPLE"), 3L)    # three rows, one dead
})

test_that("resolution keeps the row carrying the most calls", {
  cfs <- syn_cfs()
  out <- syn_run()
  markers <- setdiff(colnames(out), "SAMPLE_ID")
  raw_markers <- setdiff(colnames(cfs), SYN_META)

  for (id in c("SYN_DUP_HIGHLOW", "SYN_DUP_DISCORD", "SYN_DUP_TRIPLE")) {
    rws <- which(cfs[[SYN_ID_COL]] == id)
    counts <- vapply(rws, function(i) sum(.is_called(unlist(cfs[i, raw_markers]))),
                     integer(1))
    kept <- out[out$SAMPLE_ID == id, , drop = FALSE]
    expect_equal(nrow(kept), 1L, info = id)
    expect_equal(sum(!is.na(unlist(kept[, markers]))), max(counts),
                 info = paste("wrong row kept for", id))
  }
})

test_that("a genuinely blank cell does not outweigh real calls when picking a winner", {
  # Uncalled genotypes are spelled three ways in real exports (real NA, the
  # literal text "NA", and a blank field) and .is_called() exists precisely to
  # treat them alike. Build a pair where the "fuller-looking" row is mostly
  # blank fields and the truly fuller row spells its blanks as "NA" text, and
  # confirm the row with more real calls wins regardless of spelling.
  cfs <- syn_cfs()
  raw_markers <- setdiff(colnames(cfs), SYN_META)
  template <- cfs[cfs[[SYN_ID_COL]] == "SYN0002", ]  # a fully-called row

  row_blank <- template
  row_full <- template
  row_blank[[SYN_ID_COL]] <- "SYN_DUP_BLANKBUG"
  row_full[[SYN_ID_COL]] <- "SYN_DUP_BLANKBUG"

  # 2 real calls, the rest genuinely blank fields.
  row_blank[, setdiff(raw_markers, raw_markers[1:2])] <- ""
  # 30 real calls, the rest spelled with the literal "NA" text.
  row_full[, setdiff(raw_markers, raw_markers[1:30])] <- "NA"

  cfs2 <- rbind(cfs, row_blank, row_full)
  out <- conform_genotypes(temp_csv(cfs2), SYN_DWR, SYN_VCFS, verbose = FALSE)

  kept <- out[out$SAMPLE_ID == "SYN_DUP_BLANKBUG", , drop = FALSE]
  expect_equal(nrow(kept), 1L)
  markers <- setdiff(colnames(kept), "SAMPLE_ID")
  expect_gt(sum(.is_called(unlist(kept[, markers]))), 20)  # row_full won
})

test_that("an ID whose every row is empty is dropped entirely", {
  out <- syn_run()
  expect_false("SYN_DUP_EMPTY" %in% out$SAMPLE_ID)

  # A sample that failed but has a unique ID is NOT dropped: it is reported as
  # having no calls and passed through, because nothing else claims that ID.
  expect_true("SYN_NOCALL" %in% out$SAMPLE_ID)
})

test_that("the output has unique, non-blank SAMPLE_IDs", {
  out <- syn_run()
  expect_false(any(duplicated(out$SAMPLE_ID)))
  expect_true(all(.is_called(out$SAMPLE_ID)))
  expect_equal(nrow(out), nrow(syn_cfs()) - 6L)  # 5 dropped re-runs + 1 dead ID
})

test_that("only the genuinely discordant re-run is raised as an anomaly", {
  discordant <- grep("re-runs disagree", syn_anomalies(), value = TRUE)

  expect_length(discordant, 1L)
  expect_match(discordant, "SYN_DUP_DISCORD")
  expect_match(discordant, "\\d+ of \\d+ shared calls")
  # Re-runs of one fish that merely differ in call rate must not be reported as
  # disagreeing -- that would bury the pair that actually does.
  expect_false(any(grepl("SYN_DUP_HIGHLOW|SYN_DUP_TRIPLE", syn_anomalies())))
})

test_that("the discordance threshold is what decides", {
  # Above the fixture's ~12% discordance nothing is raised; below it, the pair
  # still is. This is the knob, not a hard-coded judgement about re-runs.
  expect_length(
    grep("re-runs disagree",
         syn_anomalies(max_replicate_discordance = 0.5)),
    0L
  )
  expect_length(
    grep("re-runs disagree",
         syn_anomalies(max_replicate_discordance = 0.05)),
    1L
  )
})

test_that("resolution reports every decision it makes", {
  txt <- syn_log()
  expect_match(txt, "SAMPLE_IDs appear more than once")
  expect_match(txt, "kept row")
  expect_match(txt, "rows removed")
  expect_match(txt, "DROPPED entirely")
})

test_that("resolution can be turned off", {
  out <- syn_run(resolve_duplicate_ids = FALSE)
  expect_equal(nrow(out), nrow(syn_cfs()))
  expect_true(any(duplicated(out$SAMPLE_ID)))
  expect_true(any(grepl("duplicate SAMPLE_IDs left in place",
                        conform_report(out)$anomalies)))

  # Duplicates are exactly why the result is not a geno_data: the rest of the
  # package assumes one row per sample and rightly refuses this table.
  expect_error(as_geno_data(out, id_col = "SAMPLE_ID", group_col = NULL),
               "Duplicate sample IDs")
})

test_that("poor and failed samples are called out by call rate", {
  txt <- syn_log()
  expect_match(txt, "samples have a call rate below")
  expect_match(txt, "samples have no calls at all")
})

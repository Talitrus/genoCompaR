synthetic_geno_df <- function() {
  data.frame(
    group     = c("G1",  "G1",  "G1",  "G2",         "G2",    "G2"),
    indiv.ID  = c("mom",  "dad",  "kid",  "unrelated1", "dup_a", "dup_b"),
    locus1.1  = c("A",   "G",   "A",   "G",           "A",     "A"),
    locus1.2  = c("A",   "G",   "G",   "G",           "G",     "G"),
    locus2.1  = c("A",   "A",   "G",   "A",           "G",     "G"),
    locus2.2  = c("G",   "G",   "G",   "A",           "G",     "G"),
    locus3.1  = c("A",   "A",   "A",   "G",           "A",     "A"),
    locus3.2  = c("A",   "A",   "A",   "G",           "A",     "A"),
    locus4.1  = c(NA,    "A",   "A",   "G",           "A",     "A"),
    locus4.2  = c(NA,    "G",   "A",   "G",           "A",     "A"),
    locus5.1  = c("A",   "G",   "C",   "A",           "C",     "T"),
    locus5.2  = c("A",   "G",   "C",   "A",           "C",     "T"),
    stringsAsFactors = FALSE
  )
}

synthetic_geno_csv <- function(df = synthetic_geno_df()) {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(df, path, row.names = FALSE, na = "NA")
  path
}

synthetic_geno <- function() {
  read_genotypes(synthetic_geno_csv())
}

example_data_path <- function() {
  system.file(
    "extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv",
    package = "genoCompaR"
  )
}

.parents_file_path <- function() {
  system.file(
    "extdata", "Parents_FRH2024_broodstock_reference_genotypes.csv",
    package = "genoCompaR"
  )
}

.unknowns_file_path <- function() {
  system.file(
    "extdata", "Unknowns_CFS_WY26_salvage_v3_03262026.csv",
    package = "genoCompaR"
  )
}

# Builds a small 3-row geno_data (offspring + two known true parents) from
# the Parents (broodstock reference) and Unknowns (salvage) files, restricted
# to their shared loci. Reads the raw CSVs directly rather than going through
# read_genotypes(), because the Unknowns file legitimately contains a handful
# of duplicate SAMPLE_IDs elsewhere (re-run samples) that read_genotypes()
# rightly rejects for whole-file loads; here we only need three specific,
# individually-unique rows. Also deliberately avoids any all-pairs scan
# across either source file (~3,200 and ~1,200 samples respectively).
.known_trio_geno <- function(offspring_id = "C260383CVP",
                              parent_ids = c("C240716FRH", "C240722FRH")) {
  p_path <- .parents_file_path()
  u_path <- .unknowns_file_path()

  read_raw <- function(path) {
    utils::read.csv(
      path,
      colClasses = "character",
      na.strings = c("NA", ""),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  parents <- read_raw(p_path)
  unknowns <- read_raw(u_path)

  p_loci <- unique(sub("\\.(1|2)$", "", setdiff(names(parents), "SAMPLE_ID")))
  u_loci <- unique(sub("\\.(1|2)$", "", setdiff(names(unknowns), "SAMPLE_ID")))
  shared_loci <- intersect(p_loci, u_loci)
  locus_cols <- as.vector(rbind(paste0(shared_loci, ".1"), paste0(shared_loci, ".2")))

  p_sub <- parents[parents$SAMPLE_ID %in% parent_ids, c("SAMPLE_ID", locus_cols)]
  u_sub <- unknowns[unknowns$SAMPLE_ID %in% offspring_id, c("SAMPLE_ID", locus_cols)]
  combined <- rbind(u_sub, p_sub)

  structure(
    tibble::as_tibble(combined),
    class = c("geno_data", class(tibble::as_tibble(combined))),
    id_col = "SAMPLE_ID",
    group_col = NA_character_,
    loci = shared_loci
  )
}

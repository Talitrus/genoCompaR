# genoCompaR

Utilities for reading and comparing wide-format Chinook salmon genotype /
haplotype call sheets — the kind produced by GT-seq panels, with one row per
sample and two allele columns per locus (`<locus>.1`, `<locus>.2`).

## Install

```r
remotes::install_github("Talitrus/genoCompaR")
```

## Quick start

```r
library(genoCompaR)

path <- system.file("extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv", package = "genoCompaR")
geno <- read_genotypes(path)
geno
#> <geno_data> 288 samples x 193 loci
```

### The `geno_data` object

`read_genotypes()` and `as_geno_data()` both return a `geno_data` object: an
ordinary tibble (subset, filter, and pipe it like any other) that also
carries three attributes the rest of the package reads instead of
hardcoding column names:

- the sample ID column (`id_col`, default `"indiv.ID"`)
- an optional grouping column (`group_col`, e.g. population/run/study group
  — default `"group"`; pass `group_col = NULL` if there is no such column)
- the locus panel (`get_loci(geno)`, without the `.1`/`.2` suffixes)

Anything that isn't a fresh CSV read — a subset of rows, a hand-picked
locus panel, a table assembled from another source — goes through
`as_geno_data()` to (re)attach these before the comparison functions will
accept it:

```r
small <- as_geno_data(as.data.frame(geno)[1:3, ], loci = get_loci(geno)[1:20])
```

## Functions

### `compare_genotypes()` — compare two samples

Loci shared, match rate, and a **bi-directional Mendelian incompatibility**
check (on by default) for parentage screening: a locus is incompatible when
the two genotypes share no allele at all, and the detail table shows which
allele(s) are missing from each side.

```r
compare_genotypes(geno, "CH-00001_S1", "CH-00002_S2")
#> Loci compared: 179   Matching: 73   Percent match: 40.8%
#> Mendelian incompatibilities: 27 (15.1% of compared loci)

# tibble output for looping over many pairs, or a plain identity check
compare_genotypes(geno, "CH-00001_S1", "CH-00002_S2", output = "tibble")
compare_genotypes(geno, "CH-00001_S1", "CH-00002_S2", mendelian = FALSE)
```

### `check_trio()` — parentage test against two known parents

Tests whether an offspring's genotype is explainable by inheriting one
allele from `parent1` and one from `parent2` at every locus. Stronger than
two separate `compare_genotypes()` duo checks, which can each pass even
when no valid one-allele-from-each assignment exists across both parents.

```r
check_trio(geno, "CH-00001_S1", "CH-00002_S2", "CH-00003_S3")
#> Compatible loci: 128   Incompatible loci: 50   Percent compatible: 71.9%
```

When the offspring and candidate parents live in separate files with
different locus panels, `merge_genotypes()` combines them on shared loci:

```r
parents <- read_genotypes(p_path, id_col = "SAMPLE_ID", group_col = NULL)
trio <- merge_genotypes(offspring, parents[parents$SAMPLE_ID %in% c("C240716FRH", "C240722FRH"), ])
check_trio(trio, "C260383CVP", "C240716FRH", "C240722FRH")
#> Compatible loci: 176   Incompatible loci: 2   Percent compatible: 98.9%
```

### `find_duplicate_genotypes()` — scan for matching/duplicate samples

Flags every pair matching at or above a threshold of their shared loci —
useful for catching re-sampled individuals or accidental re-runs.

```r
find_duplicate_genotypes(geno, threshold = 0.95, min_loci_compared = 50)
#> <geno_duplicates> 27 flagged pair(s) among 288 samples (193 loci, threshold = 95%)
```

### `genotype_match_matrix()` — full pairwise match-rate matrix

The full sample x sample match matrix (not just threshold-passing pairs),
for clustering or a heatmap.

```r
mm <- genotype_match_matrix(geno)

# Pairs sharing too few typed loci come back NA, which hclust() rejects;
# treating them as maximally distant is the conservative reading.
d <- stats::as.dist(1 - mm / 100)
d[is.na(d)] <- 1
plot(stats::hclust(d), labels = FALSE)
```

### Missingness

```r
locus_missingness(geno)   # one row per locus
sample_missingness(geno)  # one row per sample
```

### `paralog_scan()` — find collapsed paralogous loci

Flags loci called heterozygous more often than a certain rate (default: 0.25) over what is expected in Hardy-Weinberg equilibrium.

```r
ref <- read_genotypes(p_path, id_col = "SAMPLE_ID", group_col = NULL)
paralog_scan(ref)
#> <geno_paralog_scan> 193 of 221 loci scanned, 2 flagged (excess > 0.25, min_n = 50)
#> 28 locus/loci skipped: fewer than 50 complete genotypes.
#>   locus                             n alleles obs_het exp_het excess flagged
#>   NC_037104.1:56552952-56553042  1126       2   0.999   0.500  0.499 TRUE
#>   NC_037130.1:1062935-1063235    1062       2   0.930   0.498  0.432 TRUE
#>   Ots_111681-657                 3112       2   0.528   0.497  0.032 FALSE
```

### `conform_genotypes()` — convert a CFS export into the DWR layout

Reads a CFS-format genotype export, harmonizes its marker names, validates
every allele against reference VCFs, resolves repeated sample IDs, and returns
a table matching a DWR reference file's exact column layout — the format
CKMRsim and CKMR-shiny expect.

```r
out <- conform_genotypes(
  "WY26_salvage_combined_genotypes.csv",       # CFS export
  "FRH2025_broodstock_genos.csv",              # DWR reference (column template)
  system.file("extdata", "reference_VCFs",
              c("full.vcf", "Chinook_FullPanel_VGLL3Six6LFARWRAP.vcf"),
              package = "genoCompaR")
)
#> ... offset report, duplicate-ID report, call rates, coverage, anomalies

conform_report(out)$anomalies   # what needs a human decision
conform_report(out)$offsets     # per-locus shift summary
conform_report(out)$het_excess  # the paralog_scan() result

write.csv(out, "conformed.csv", row.names = FALSE, na = "")
```

Two behaviors are **off by default** because they discard or rewrite data:
`repair_offsets`, which moves genotypes the scan diagnosed as shifted back to
the loci they belong to, and `null_invalid_alleles`, which blanks alleles the
VCFs say are impossible. Run once with both off, read the report, and turn one
on only once you accept its diagnosis.

The result carries duplicate sample IDs when `resolve_duplicate_ids = FALSE`,
so it is a plain table rather than a `geno_data`. Convert it to compare
samples:

```r
geno <- as_geno_data(out, id_col = "SAMPLE_ID", group_col = NULL)
```

## Function reference

| Function                     | Purpose                                                            |
|-------------------------------|--------------------------------------------------------------------|
| `read_genotypes()`            | Load a wide genotype/haplotype CSV into a `geno_data` object       |
| `as_geno_data()`              | Turn an existing data frame / subset into a `geno_data` object     |
| `merge_genotypes()`           | Combine two genotype sets on their shared loci                     |
| `compare_genotypes()`         | Pairwise duo comparison + optional bi-directional Mendelian check  |
| `check_trio()`                | Trio-based Mendelian parentage test (offspring + two parents)      |
| `locus_detail()`              | Full per-locus table from a `compare_genotypes()`/`check_trio()` result |
| `find_duplicate_genotypes()`  | Dataset-wide scan for matching / duplicate genotypes               |
| `genotype_match_matrix()`     | Full sample x sample match-rate matrix, for clustering/heatmaps    |
| `locus_missingness()`         | Per-locus missing-call summary                                     |
| `sample_missingness()`        | Per-sample missing-call summary                                    |
| `paralog_scan()`              | Heterozygote-excess scan for collapsed paralogous loci             |
| `conform_genotypes()`         | Convert a CFS export into the DWR haplotype layout, VCF-validated  |
| `conform_report()`            | QC report from a `conform_genotypes()` run                         |
| `get_loci()`                  | Locus names carried by a `geno_data` object                        |

## Example data

Real Chinook haplotype call sheets in `inst/extdata/`:

- `AFTC_ILC_P134_full_filtered_haplotypes.csv` — 288 samples, 193 loci.
- `Parents_FRH2024_broodstock_reference_genotypes.csv` — 3,205 broodstock
  reference samples, 221 loci. Use `id_col = "SAMPLE_ID"`.
- `Unknowns_CFS_WY26_salvage_v3_03262026.csv` — 1,236 field-salvage
  samples, 193 loci, `id_col = "SAMPLE_ID"`. Contains a few duplicate
  `SAMPLE_ID`s (re-run samples), so a whole-file `read_genotypes()` call
  errors by design — read it with `utils::read.csv()`, subset to the rows
  you need, and pass those to `as_geno_data()`.

The Parents/Unknowns files contain one known true trio (see
`tests/testthat/test-known-trio.R`): offspring `C260383CVP` (Unknowns) has
parents `C240716FRH` and `C240722FRH` (Parents). The Unknowns file is itself
the output of `conform_genotypes()`.

Also shipped:

- `inst/extdata/reference_VCFs/` — the two panel VCFs `conform_genotypes()`
  validates alleles against (203 loci between them).
- `inst/extdata/synthetic/` — five small invented files, ~50 kB, that the
  whole conversion test suite runs on; see that directory's `README.md` for
  what each fault in them is there to test, and
  `data-raw/make_synthetic_fixtures.R` to regenerate them.

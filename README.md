# genoCompaR

Utilities for reading and comparing wide-format Chinook salmon genotype /
haplotype call sheets — the kind produced by GT-seq panels, with one row per
sample and two allele columns per locus (`<locus>.1`, `<locus>.2`).

## Install

```r
devtools::install(".")
```

## Quick start

```r
library(genoCompaR)

path <- system.file("extdata", "AFTC_ILC_P134_full_filtered_haplotypes.csv", package = "genoCompaR")
geno <- read_genotypes(path)
geno
#> <geno_data> 288 samples x 193 loci
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
parents `C240716FRH` and `C240722FRH` (Parents).

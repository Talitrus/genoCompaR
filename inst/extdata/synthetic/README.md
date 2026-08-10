# Synthetic fixtures

Five files, ~50 kB, no real fish. They reproduce in miniature every fault and
every invariant the `conform_genotypes()` test suite asserts, so the suite runs
with no CDFW/DWR genotypes present at all. Reach them with

```r
system.file("extdata", "synthetic", "synthetic_cfs_clean.csv",
            package = "genoCompaR")
```

or, in the test suite, through the `SYN_*` constants in
`tests/testthat/helper-conform.R`.

Regenerate with:

```
Rscript data-raw/make_synthetic_fixtures.R
```

The generator is seeded, so output is stable run to run. Changing it changes the
fixtures; `tests/testthat/test-synthetic-fixtures.R` is what proves a
regenerated set still carries every property below — if an edit quietly stops
the data from containing one of these faults, a test goes red rather than the
suite going green on data that no longer tests anything.

## Files

| file | what it is |
|---|---|
| `synthetic_full.vcf` | 39 SNP records over 18 loci, plus 4 indel records (2 inside a SNP locus, 2 forming a 19th locus that has to vanish) |
| `synthetic_panel.vcf` | 3 loci, 2 of them repeated from `synthetic_full.vcf` at identical `(CHROM, POS)` |
| `synthetic_dwr_reference.csv` | 80 reference samples × 22 loci, DWR naming, columns sorted by locus |
| `synthetic_cfs_clean.csv` | 43 rows / 38 IDs × 20 loci, CFS naming, correct header order |
| `synthetic_cfs_offset.csv` | byte-identical data, header carrying the 20260313 fault |

The panel spans haplotype widths 1–5 and both marker-name conventions. 19 loci
are described by the VCFs, 20 are on the CFS panel (`tag_id_999` is on neither
VCF), and 2 are DWR-only.

## What each feature is there to test

### Marker naming — `transform_name()`

| locus (CFS → DWR) | rule |
|---|---|
| `tag_id_101` → `tag_id_101.1` | allele suffix shifts; name unchanged |
| `NC_037099.1.62937268.62937373` → `NC_037099.1:62937268-62937373` | genomic coordinates; the `.1` assembly version must survive |
| `Ots_111312.435` → `Ots_111312-435` | dots become dashes |
| `OkiOts_120255.113` → `OkiOts_120255-113` | same, `OkiOts_` prefix |
| `Ots_BMP.2.SNP1` → `Ots_BMP-2-SNP1` | **every** dot, not just the first — the `CFS_to_DWR2.R` regression that left the locus silently blank |
| `Ots_u07.07.161` → `Ots_u07-07.161` | the exception: only the first dot, so the DWR name keeps a literal dot |

### VCF loading

- `NC_037099.1:62937268-62937373` and `NC_037130.1:1062935-1063235` appear in
  **both** VCFs at identical `(CHROM, POS)`. Without dedup their widths double
  from 2 to 4 and every real allele at them turns invalid. (In `full.vcf` this
  is 8 repeated `NC_` loci.)
- `NC_037124.1:12281207-12281551` has 2 SNP sites **and** 2 indel records.
  Counting the indels makes it 4 wide and admits multi-base "alleles".
- `NC_037124.1:12270118-12270418` is nothing but indel records, so after
  filtering it has to disappear from the panel rather than appear zero-wide.
- `tag_id_104` site 2 is `REF A, ALT G,AT`. The site stays; `AT` must not
  become a permitted base.

### Reference consistency

- Every DWR locus width equals its VCF site count, so when CFS disagrees with
  them it is CFS that is wrong.
- `tag_id_106` carries the allele **`TCAGA`**, present in 5 reference samples
  and in 2 CFS samples, and absent from every VCF ALT list — its first base is
  `T` where the VCF permits only `A`/`G`. This is the single-sample-VCF problem:
  `build_permitted()` has to widen the permitted set with what the reference
  actually shows, or genuine variation reads as corruption and
  `NULL_INVALID_ALLELES` deletes it. It stays well under 1% of reference
  alleles, so a strict-whitelist regression is still detectable as a small
  non-zero.

### The header-order fault

`synthetic_cfs_offset.csv` holds exactly the same genotype matrix as
`synthetic_cfs_clean.csv`; only the header line is permuted. `tag_id_101` is
listed at the **end** of the `tag_id` block instead of the front — the shape of
the 20260313 export, where `tag_id_5720` did the same and mislabelled 18 loci.

Reading the result: the value under label `tag_id_102` is really `tag_id_101`'s
genotype, so every genotype from there on belongs one locus to the left. The
scan reports it as **shift −1 across 9 consecutive loci in ~81% of samples**,
with one locus (`tag_id_102`, holding data from 9 loci away) left as
unexplained-in-place.

Two properties make this a real test rather than a width check:

- Widths differ at **every adjacent pair** in the block (2, 3, 1, 4, 2, 5, 2, 3,
  1, 2), so most of the slip is caught outright.
- `tag_id_110` and `tag_id_101` are **both 2 wide**. The slip between them is
  invisible to a length test and is caught only because their permitted bases at
  site 1 are disjoint (`{C,T}` vs `{A,G}`). This is the case that motivated
  validating alleles against the VCF instead of measuring lengths.

Measured behaviour, default configuration:

| | `REPAIR_OFFSETS = FALSE` | `TRUE` |
|---|---|---|
| genotypes retained | 744 | 1246 |
| alleles still invalid | 70 | 14 |

Some invalid alleles survive repair by design: the genotype shifted off the
front of the block is unrecoverable, and `NULL_INVALID_ALLELES` is off, so
nothing is deleted silently. Setting it to `TRUE` takes the count to 0.

### Duplicate `SAMPLE_ID`s

All four shapes are present. Rows sharing an ID are re-runs of one fish, so they
start from one genotype and differ only in how much of it was called — drawing
them independently would make every duplicate look like a replicate conflict and
bury the pair that genuinely is one.

| ID | rows | expected outcome |
|---|---|---|
| `SYN_DUP_HIGHLOW` | 2 | failed run + good run; keep the good one, no anomaly |
| `SYN_DUP_EMPTY` | 2 | both empty; **ID dropped entirely** |
| `SYN_DUP_DISCORD` | 2 | two good runs disagreeing at ~11% of shared calls; keep the fuller row **and** raise an anomaly |
| `SYN_DUP_TRIPLE` | 3 | one dead + two live; keep the fullest, no anomaly |

`SYN_DUP_DISCORD` is the only anomaly a clean run reports.

### Call rates and coverage

| sample / locus | why it is there |
|---|---|
| `SYN_LOWCALL` | ~35% call rate → the `LOW_CALL_RATE` warning |
| `SYN_NOCALL` | no calls, but a **unique** ID → reported as empty and passed through, not dropped |
| `tag_id_999` | in both CSVs, in neither VCF → width from the reference only, alleles unchecked; must appear in the COVERAGE section |
| `winter-greb1l-diagnostic_09_*`, `_10_*` | DWR-only → present in the output and entirely empty |

The CFS file also contains both spellings of an uncalled genotype — the literal
string `NA` and a genuinely empty cell — since exports have used both and
`is_called()` has to treat them alike.

### Where the measured numbers are asserted

The counts above are not decoration — they are checked. `744`/`70` and
`1246`/`14` are asserted in `tests/testthat/test-conform-offsets.R`, the
duplicate-ID outcomes in `test-conform-duplicates.R`, the VCF widths and dedup
in `test-vcf-loading.R`, and the `tag_id_106` widening gap in
`test-allele-validation.R`. If an edit to the generator quietly stops the data
from carrying one of these properties, a test goes red rather than the suite
going green on data that no longer tests anything.

One fault is deliberately **not** here: there is no collapsed-paralog locus in
this set. `paralog_scan()` is tested against the two real ones in
`inst/extdata/Parents_FRH2024_broodstock_reference_genotypes.csv`, and against
an injected locus for the `conform_genotypes()` reporting path — see
`tests/testthat/test-paralog_scan.R`.

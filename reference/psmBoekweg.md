# Example peptide-spectrum match data

`psmBoekweg` is a `PSM` object containing 10,000 peptide-spectrum
matches randomly sampled from a larger dataset. The PSMs were identified
using the Sage search engine.

## Format

A `PSM` object with 10,000 rows. Each PSM contains the following
metadata columns (among others):

- `psm_id`:

  Unique identifier for the PSM

- `peptide`:

  Peptide sequence

- `proteins`:

  Protein identifiers

- `filename`:

  Source mzML file name

- `scannr`:

  Scan number

- `rank`:

  Rank of the PSM

- `charge`:

  Precursor charge state

- `hyperscore`:

  Sage hyperscore

- `peptide_q`:

  Peptide-level q-value (FDR)

- `protein_q`:

  Protein-level q-value (FDR)

- `rt`:

  Retention time (seconds)

## Details

This is an example subset of the data available through
`?MsDataHub::Boekweg2022()`.

## See also

- [spBoekweg](https://rformassspectrometry.github.io/PSMatch/reference/spBoekweg.md)
  for the corresponding `Spectra` object.

- [`PSM()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md)
  for the PSM class constructor.

- The `inst/scripts/make_BoekwegData.R` script documenting data
  creation.

## Examples

``` r
data(psmBoekweg)
psmBoekweg
#> PSM with 700 rows and 40 columns.
#> Spectra: 700 unique
#>   db: 572 target, 128 decoy
#>   ranks: 1:700 
#> names(40): psm_id peptide ... protein_q ms2_intensity

## Number of PSMs
length(psmBoekweg)
#> [1] 40

## Access peptide sequences
head(psmBoekweg$peptide)
#> [1] "VATAQDDITGDGTTSNVLIIGELLK"   "ASASGSGAQVGGPISSGSSASSVTVTR"
#> [3] "DLGLSESGEDVNAAILDESGK"       "FDGALNVDLTEFQTNLVPYPR"      
#> [5] "MSVQPTVSLGGFEITPPVVLR"       "LGQSDPAPLQHQMDIYQK"         

## Filter to rank 1 PSMs only
psmRank1 <- filterPsmRank(psmBoekweg)
#> Removed 0 PSMs with rank > 1.
length(psmRank1)
#> [1] 40

## Summary of q-values
summary(psmBoekweg$peptide_q)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 4.338e-05 4.338e-05 4.338e-05 7.034e-02 1.200e-01 3.598e-01 
```

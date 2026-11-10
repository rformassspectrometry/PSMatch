# Filter out unreliable PSMs.

Functions to filter out PSMs matching. The PSMs should be stored in a
`PSM` such as those produced by
[`PSM()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).

- `filterPsmDecoy()` filters out decoy PSMs, i.e. those annotated as
  `isDecoy`.

&nbsp;

- `filterPsmRank()` filters out PSMs of rank \> 1.

&nbsp;

- `filterPsmShared()` filters out shared PSMs, i.e. those that match
  multiple proteins.

&nbsp;

- `filterPsmFdr()` filters out PSMs based on their FDR.

## Usage

``` r
filterPSMs(
  x,
  decoy = psmVariables(x)["decoy"],
  rank = psmVariables(x)["rank"],
  protein = psmVariables(x)["protein"],
  spectrum = psmVariables(x)["spectrum"],
  peptide = psmVariables(x)["peptide"],
  verbose = TRUE
)

filterPsmDecoy(x, decoy = psmVariables(x)["decoy"], verbose = TRUE)

filterPsmRank(x, rank = psmVariables(x)["rank"], verbose = TRUE)

filterPsmShared(
  x,
  protein = psmVariables(x)["protein"],
  peptide = psmVariables(x)["peptide"],
  verbose = TRUE
)

filterPsmFdr(x, FDR = 0.05, fdr = psmVariables(x)["fdr"], verbose = TRUE)
```

## Arguments

- x:

  An instance of class `PSM`.

- decoy:

  `character(1)` with the column name specifying whether entries match
  the decoy database or not. Default is the `decoy` PSM variable as
  defined in
  [`psmVariables()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).
  The column should be a `logical` and only PSMs holding a `FALSE` are
  retained. Filtering is ignored if set to `NULL` or `NA`.

- rank:

  `character(1)` with the column name holding the rank of the PSM.
  Default is the `rank` PSM variable as defined in
  [`psmVariables()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).
  This column should be a `numeric` and only PSMs having rank equal to 1
  are retained. Filtering is ignored if set to `NULL` or `NA`.

- protein:

  `character(1)` with the column name holding the protein (groups)
  protein. Default is the `protein` PSM variable as defined in
  [`psmVariables()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).
  Filtering is ignored if set to `NULL` or `NA`.

- spectrum:

  `character(1)` with the name of the spectrum identifier column.
  Default is the `spectrum` PSM variable as defined in
  [`psmVariables()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).
  Filtering is ignored if set to `NULL` or `NA`.

- peptide:

  `character(1)` with the name of the peptide identifier column. Default
  is the `peptide` PSM variable as defined in
  [`psmVariables()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).
  Filtering is ignored if set to `NULL` or `NA`.

- verbose:

  `logical(1)` setting the verbosity flag.

- FDR:

  `numeric(1)` to be used to filter based on the `fdr` variable. Default
  is 0.05.

- fdr:

  `character(1)` variable name that defines that defines the spectrum
  FDR (or any similar/relevant metric that can be used for filtering).
  This value isn't set by default as it depends on the search engine and
  application. Default is `NA`.

## Value

A new filtered `PSM` object with the same columns as the input `x`.

## Author

Laurent Gatto

## Examples

``` r
f <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache

id <- PSM(f)
filterPSMs(id)
#> Starting with 5802 PSMs:
#> Removed 2896 decoy hits.
#> Removed 155 PSMs with rank > 1.
#> Removed 85 shared peptides.
#> 2666 PSMs left.
#> PSM with 2666 rows and 35 columns.
#> names(35): sequence spectrumID ... subReplacementResidue subLocation
```

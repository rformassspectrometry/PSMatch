# Describe protein and peptide compositions

It is important to explore PSM results prior to any further downstream
analysies. Two functions, that work on
[`PSM()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md)
and
[`ConnectedComponents()`](https://rformassspectrometry.github.io/PSMatch/reference/ConnectedComponents.md)
objects can be used for this:

- The `describeProteins()` function describe protein composition in
  terms of unique and shared peptides.

- The `describePeptides()` function describe unique/shared peptide
  composition.

## Usage

``` r
describeProteins(object, ...)

describePeptides(object, ...)
```

## Arguments

- object:

  Either an instance of class `Matrix`,
  [`PSM()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md)
  or
  [`ConnectedComponents()`](https://rformassspectrometry.github.io/PSMatch/reference/ConnectedComponents.md).

- ...:

  Additional arguments passed to
  [`makeAdjacencyMatrix()`](https://rformassspectrometry.github.io/PSMatch/reference/adjacencyMatrix.md).

## Value

`describePeptides()` invisibly return the table of unique and shared
peptides. `describeProteins()` invisibly returns a `data.frame` with
logicals indicating the unique/shared peptide composition of proteins.
Both functions are used for their side effects of printing a short
descriptive output about peptides and proteins.

## Examples

``` r
f <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache

psm <- PSM(f) |>
       filterPsmDecoy() |>
       filterPsmRank()
#> Removed 2896 decoy hits.
#> Removed 155 PSMs with rank > 1.

describePeptides(psm)
#> 2357 peptides composed of
#>  unique peptides: 2324
#>  shared peptides (among protein):
#>   28(2) 3(3) 2(4)
describeProteins(psm)
#> 1504 proteins composed of
#>  only unique peptides: 1455
#>  only shared peptides: 38
#>  unique and shared peptides: 11
```

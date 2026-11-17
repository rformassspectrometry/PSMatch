# findBetterMatch

Try to find a better peptide-spectrum match than the current PSM.

## Usage

``` r
findBetterMatch(
  x,
  peptideVariable = "peptide",
  modificationsList = NULL,
  maxMods = 2,
  improvedBy = c("maxScore", "maxCount"),
  convertToStyle = "name"
)
```

## Arguments

- x:

  `Spectra()` object containing identifications.

- peptideVariable:

  `character(1L)` pointing to the peptide sequence variable in `x`.

- modificationsList:

  [`character()`](https://rdrr.io/r/base/character.html) with all
  modifications to try out. The modifications need to be present in the
  [PTMods::modifications](https://rdrr.io/pkg/PTMods/man/modifications.html)
  dataset from the unimod database. Note that the more modifications are
  included, the more combinations are possible and the longer it takes
  to run.

- maxMods:

  `integer(1L)` indicating the maximum modifications allowed at the same
  time on one sequence.

- improvedBy:

  `character(1L)`: either `"maxScore"` or `"maxCount"`, defining how to
  choose the better match. `"maxScore"` keeps the match with the highest
  `ndotproduct` and `"maxCount"` keeps the one with the most matched
  peaks.

- convertToStyle:

  `character(1L)` to choose from those available in
  [`PTMods::convertAnnotation()`](https://rdrr.io/pkg/PTMods/man/convertAnnotation.html).
  By default, the unimod "name" is given.

## Value

A named [`character()`](https://rdrr.io/r/base/character.html) vector of
length equal to the number of spectra with identifications. The values
are the best matches found using the `modificationsList` and the names
are the indices for those sequences in the innitial Spectra object `x`.

## Details

This function is meant to help finding a better peptide-spectrum match
for identifications that fail certain checks (i.e.
[`checkOverlap()`](https://rformassspectrometry.github.io/PSMatch/reference/validatePSM.md)
has returned false) or if the PSM's FDR doesn't qualify it as a
confident identification. A peptide sequence is returned (with or
without modifications) for which the greatest number of matched peaks/a
better ndotproduct (see
[`MsCoreUtils::ndotproduct()`](https://rdrr.io/pkg/MsCoreUtils/man/distance.html)
for more details) was obtained.

## Author

Guillaume Deflandre

## Examples

``` r
if (FALSE) { # \dontrun{
## Load the data: a Spectra object with identifications
data(spBoekweg)
data(psmBoekweg)
head(psmBoekweg$pkey <- paste0(basename(psmBoekweg$filename),
                               sub("^.+scan=", "::", psmBoekweg$scannr)))
head(spBoekweg$pkey <- paste0(basename(spBoekweg$dataOrigin),
                              sub("^.+scan=", "::", spBoekweg$spectrumId)))
sp <- Spectra::joinSpectraData(spBoekweg, psmBoekweg, by.x = "pkey")

## Run with stock modifications (see PSMatch:::.ptmNames)
findBetterMatch(sp[30], "peptide", maxMods = 1)

## Run with custom modifications (need to be present in
## `modifications$Name`). Notice that for this example, the identified match
## includes an oxidation. We are removing it to show that findBetterMatch()
## also identifies the same sequence.
sp[30]$sequence <- getCanonicalSequence(sp[30]$peptide)
findBetterMatch(sp[30], "sequence", maxMods = 1)

## Get the best match based on the number of matched peaks
findBetterMatch(sp[30], "sequence", maxMods = 1, improvedBy = "maxCount")
} # }
```

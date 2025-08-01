# PSMatch 1.13

## PSMatch 1.13.3

- Nothing yet

## PSMatch 1.13.2

- Add CITATION to the pre-print (https://doi.org/10.31219/osf.io/62v9p_v2)
- Add possibility to split protein groups in PSM data with `makeAdjacencyMatrix`
(see [issue #30](https://github.com/rformassspectrometry/PSMatch/issues/30))
- Removed `showDetails()` from `setMethod("show", "PSM")` as discussed in
[issue #30](https://github.com/rformassspectrometry/PSMatch/issues/30)
- Add `USI` parameter to `plotSpectraPTM`.
- Correct `filterPsmFdr` output message
- Update Fragments vignette.
- Improve `labelFragments()`runtime
(see [issue #25](https://github.com/rformassspectrometry/PSMatch/issues/25)).
- Add Selenocysteine and Pyrrolysine to `getAminoAcids()`.

## PSMatch 1.13.1

- Refactor and improve `plotSpectraPTM()`.

## PSMatch 1.13.0

- New Bioconductor devel.

# PSMatch 1.11

## PSMatch 1.11.5

- Correct `plotSpectraPTM` annotations.

## PSMatch 1.11.4

- Added `plotSpectraPTM`: a plotting function to visualise
  post-translational modifications.

## PSMatch 1.11.3

- Deprecated `addFragments`. The use of `labelFragments` is endorsed
  instead.  See [PR #20](https://github.com/rformassspectrometry/PSMatch/pull/20).

## PSMatch 1.11.2

- Replace `calculateFragments` with `calculateFragments2`. See
  [PR #19](https://github.com/rformassspectrometry/PSMatch/pull/19).

## PSMatch 1.11.1

- New `calculateFragments2` function includes fixed and variable
  modifications to fragments ions. See [PR #16](https://github.com/rformassspectrometry/PSMatch/pull/16).

## PSMatch 1.11.0

- New devel version

# PSMatch 1.9

## PSMatch 1.9.1

- Fix check errors.

## PSMatch 1.9.0

- New Bioc devel.

# PSMatch 1.7

## PSMatch 1.7.2

- Fix connected component dim names in `show()`.

## PSMatch 1.7.1

- In `addFragments()` use `...` to pass parameters to
  `calculateFragments()`.

## PSMatch 1.7.0

- New Bioc devel.

# PSMatch 1.5

## PSMatch 1.5.0

- New Bioc devel.

# PSMatch 1.3

## PSMatch 1.3.3

- New `fdr` variable (default is always `NA_character_` for now) that
  defines the spectrum FDR (or any similar/relevant metric that can be
  used for filtering - see next item).
- New `filterPsmFdr()` function that filters based on the `fdr`
  variable.

## PSMatch 1.3.2

- Specific `Matrix::rowSums()` to fix error in example.

## PSMatch 1.3.1

- Fix type in vignette.

# PSMatch 1.0

## PSMatch 1.0.0

- First Bioconductor release.

# PSMatch 0.99

## Changes in 0.99.5

- Fix *mz* calculation in `calculateFragments` for neutral losses with
  a charge > 1 (ported from `MSnbase` - see [issue
  573](https://github.com/lgatto/MSnbase/issues/573)).

## Changes in 0.99.4

- Set seed in the ConnectedComponents unit test to stop random errors
  after clustering.

## Changes in 0.99.3

- Fix bug in `describePeptides()` (close #11).

## Changes in 0.99.2

- Describe the `ConnectedComponents()` return value.
- Add/update installation instructions.

## Changes in 0.99.1

- Fix typo and improve documentation.

## Changes in 0.99.0

- Prepare package for Bioconductor submission.
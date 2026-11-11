# Changelog

## PSMatch 1.17

### PSMatch 1.17.3

- Change license to GPL-3.

### PSMatch 1.17.2

- Added
  [`validatePSM()`](https://rformassspectrometry.github.io/PSMatch/reference/validatePSM.md)
  and its corresponding vignette and tests.
- Added data bundled with the package: `spBoekweg` and `psmBoekweg`
- Added paper reference in CITATION

### PSMatch 1.17.1

- Added parameter `z` and `allCharges` in
  [`plotSpectraPTM()`](https://rformassspectrometry.github.io/PSMatch/reference/plotSpectraPTM.md)
  and
  [`labelFragments()`](https://rformassspectrometry.github.io/PSMatch/reference/labelFragments.md)
  and set it to `1:precursorCharge(x)` by default when it used to only
  display charge 1 before. See this
  [PR](https://github.com/rformassspectrometry/PSMatch/pull/44).

### PSMatch 1.17.0

- New Bioconductor devel.

## PSMatch 1.15

### PSMatch 1.15.3

- Adjusted documentation on `addCarbamidomethyl = TRUE` in
  [`calculateFragments()`](https://rformassspectrometry.github.io/PSMatch/reference/calculateFragments.md)
  set by default.
- Corrected
  [`plotSpectraPTM()`](https://rformassspectrometry.github.io/PSMatch/reference/plotSpectraPTM.md)
  relative `PTMods` dependencies, identifications are highlighted in
  bold in the USI. Added parameters to call addFixed and addVariable
  within
  [`plotSpectraPTM()`](https://rformassspectrometry.github.io/PSMatch/reference/plotSpectraPTM.md).
- Add PTMods dependency and thus positional modifications in
  [`calculateFragments()`](https://rformassspectrometry.github.io/PSMatch/reference/calculateFragments.md)
  (see [issue
  38](https://github.com/rformassspectrometry/PSMatch/issues/38))

### PSMatch 1.15.2

- Use `MsDataHub` instead of `msdata`.

### PSMatch 1.15.1

- Fix minor check notes (Rd link and IRanges import)

### PSMatch 1.15.0

- New devel version

## PSMatch 1.13

### PSMatch 1.13.3

- Nothing yet

### PSMatch 1.13.2

- Add CITATION to the pre-print
  (<https://doi.org/10.31219/osf.io/62v9p_v2>)
- Add possibility to split protein groups in PSM data with
  `makeAdjacencyMatrix` (see
  [issue](https://github.com/rformassspectrometry/PSMatch/issues/30)
  [\#30](https://github.com/RforMassSpectrometry/PSMatch/issues/30))
- Removed `showDetails()` from `setMethod("show", "PSM")` as discussed
  in [issue](https://github.com/rformassspectrometry/PSMatch/issues/30)
  [\#30](https://github.com/RforMassSpectrometry/PSMatch/issues/30)
- Add `USI` parameter to `plotSpectraPTM`.
- Correct `filterPsmFdr` output message
- Update Fragments vignette.
- Improve
  [`labelFragments()`](https://rformassspectrometry.github.io/PSMatch/reference/labelFragments.md)runtime
  (see
  [issue](https://github.com/rformassspectrometry/PSMatch/issues/25)
  [\#25](https://github.com/RforMassSpectrometry/PSMatch/issues/25)).
- Add Selenocysteine and Pyrrolysine to
  [`getAminoAcids()`](https://rformassspectrometry.github.io/PSMatch/reference/getAminoAcids.md).

### PSMatch 1.13.1

- Refactor and improve
  [`plotSpectraPTM()`](https://rformassspectrometry.github.io/PSMatch/reference/plotSpectraPTM.md).

### PSMatch 1.13.0

- New Bioconductor devel.

## PSMatch 1.11

### PSMatch 1.11.5

- Correct `plotSpectraPTM` annotations.

### PSMatch 1.11.4

- Added `plotSpectraPTM`: a plotting function to visualise
  post-translational modifications.

### PSMatch 1.11.3

- Deprecated `addFragments`. The use of `labelFragments` is endorsed
  instead. See
  [PR](https://github.com/rformassspectrometry/PSMatch/pull/20)
  [\#20](https://github.com/RforMassSpectrometry/PSMatch/issues/20).

### PSMatch 1.11.2

- Replace `calculateFragments` with `calculateFragments2`. See
  [PR](https://github.com/rformassspectrometry/PSMatch/pull/19)
  [\#19](https://github.com/RforMassSpectrometry/PSMatch/issues/19).

### PSMatch 1.11.1

- New `calculateFragments2` function includes fixed and variable
  modifications to fragments ions. See
  [PR](https://github.com/rformassspectrometry/PSMatch/pull/16)
  [\#16](https://github.com/RforMassSpectrometry/PSMatch/issues/16).

### PSMatch 1.11.0

- New devel version

## PSMatch 1.9

### PSMatch 1.9.1

- Fix check errors.

### PSMatch 1.9.0

- New Bioc devel.

## PSMatch 1.7

### PSMatch 1.7.2

- Fix connected component dim names in
  [`show()`](https://rdrr.io/r/methods/show.html).

### PSMatch 1.7.1

- In `addFragments()` use `...` to pass parameters to
  [`calculateFragments()`](https://rformassspectrometry.github.io/PSMatch/reference/calculateFragments.md).

### PSMatch 1.7.0

- New Bioc devel.

## PSMatch 1.5

### PSMatch 1.5.0

- New Bioc devel.

## PSMatch 1.3

### PSMatch 1.3.3

- New `fdr` variable (default is always `NA_character_` for now) that
  defines the spectrum FDR (or any similar/relevant metric that can be
  used for filtering - see next item).
- New
  [`filterPsmFdr()`](https://rformassspectrometry.github.io/PSMatch/reference/filterPSMs.md)
  function that filters based on the `fdr` variable.

### PSMatch 1.3.2

- Specific
  [`Matrix::rowSums()`](https://rdrr.io/pkg/Matrix/man/colSums-methods.html)
  to fix error in example.

### PSMatch 1.3.1

- Fix type in vignette.

## PSMatch 1.0

### PSMatch 1.0.0

- First Bioconductor release.

## PSMatch 0.99

### Changes in 0.99.5

- Fix *mz* calculation in `calculateFragments` for neutral losses with a
  charge \> 1 (ported from `MSnbase` - see [issue
  573](https://github.com/lgatto/MSnbase/issues/573)).

### Changes in 0.99.4

- Set seed in the ConnectedComponents unit test to stop random errors
  after clustering.

### Changes in 0.99.3

- Fix bug in
  [`describePeptides()`](https://rformassspectrometry.github.io/PSMatch/reference/describeProteins.md)
  (close
  [\#11](https://github.com/RforMassSpectrometry/PSMatch/issues/11)).

### Changes in 0.99.2

- Describe the
  [`ConnectedComponents()`](https://rformassspectrometry.github.io/PSMatch/reference/ConnectedComponents.md)
  return value.
- Add/update installation instructions.

### Changes in 0.99.1

- Fix typo and improve documentation.

### Changes in 0.99.0

- Prepare package for Bioconductor submission.

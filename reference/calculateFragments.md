# Calculate ions produced by fragmentation with variable modifications

This method calculates a-, b-, c-, x-, y- and z-ions produced by
fragmentation.

Available methods

- The default method with signature `sequence = "character"` and
  `object = "missing"` calculates the theoretical fragments for a
  peptide sequence. It returns a `data.frame` with the columns `mz`,
  `ion`, `type`, `pos`, `z`, `seq` and `peptide`.

- Additional method can be defined that will adapt their behaviour based
  on spectra defined in `object`. See for example the MSnbase package
  that implements a method for objects of class `Spectrum2`.

## Usage

``` r
# S4 method for class 'character,missing'
calculateFragments(
  sequence,
  type = c("b", "y"),
  z = 1,
  fixed_modifications = NULL,
  variable_modifications = NULL,
  addCarbamidomethyl = TRUE,
  max_mods = Inf,
  neutralLoss = defaultNeutralLoss(),
  verbose = TRUE
)
```

## Arguments

- sequence:

  [`character()`](https://rdrr.io/r/base/character.html) providing a
  peptide sequence. If positional modifications are included in the
  sequence, variable modifications may not be used. See examples below
  for more detail.

- type:

  `character` vector of target ions; possible values:
  `c("a", "b", "c", "x", "y", "z")`. Default is `type = c("b", "y")`.

- z:

  numeric with a desired charge state; default is 1.

- fixed_modifications:

  Deprecated parameter. Please use
  [`PTMods::addFixedModifications()`](https://rdrr.io/pkg/PTMods/man/addFixedModifications.html)
  to generate sequences with positional modifications instead. Named
  `numeric` or `character`. If a `character` is given, values must be in
  UniMod name or UniMod ID format (e.g. `"Phospho"`, `"UNIMOD:21"`). The
  annotation style of the values is preserved in the output. Specifies
  which fixed modifications are applied to which amino acids.

- variable_modifications:

  Deprecated parameter. Please use
  [`PTMods::addVariableModifications()`](https://rdrr.io/pkg/PTMods/man/addVariableModifications.html)
  to generate sequences with positional modifications instead. Named
  `numeric` or `character`. If a `character` is given, values must be in
  UniMod name or UniMod ID format (e.g. `"Phospho"`, `"UNIMOD:21"`). The
  annotation style of the values is preserved in the output. Specifies
  which variable modifications are used on which amino acids.

- addCarbamidomethyl:

  logical(1L) set to `TRUE` by default. Applies carbamidomethylation as
  a fixed modification unless `fixed_modifications` is not `NULL` or if
  a carbamidomethyl is already present in the given sequences. It is
  strongly suggested to rely on
  [`PTMods::addFixedModifications()`](https://rdrr.io/pkg/PTMods/man/addFixedModifications.html)
  instead.

- max_mods:

  Deprecated parameter. Please use in combination with
  PTMods::addVariableModifications() instead. A numeric indicating the
  maximum number of variable modifications allowed on the sequence at
  once. Does not include fixed modifications. Default value is positive
  infinity.

- neutralLoss:

  `list`, it has to have two named elments, namely `water` and `ammonia`
  that contain a `character` vector which type of neutral loss should be
  calculated. Currently neutral loss on the C terminal `"Cterm"`, at the
  amino acids `c("D", "E", "S", "T")` for `"water"` (shown with an `_`)
  and `c("K", "N", "Q", "R")` for `"ammonia"` (shown with an `*`) are
  supported.

      There is a helper function `defaultNeutralLoss()` that returns
      the correct list. It has two arguments `disableWaterLoss` and
      `disableAmmoniaLoss` to remove single neutral loss options. See
      the example section for use cases.

- verbose:

  `logical(1)`. Deprecated parameter. If `TRUE` (default) the used
  modifications are printed.

## Value

A `data.frame` showing all the ions produced by fragmentation with all
possible combinations of modifications. The used variable modifications
are displayed in the `peptide` column through the use of amino acids
followed by the modification within brackets. Fixed modifications are
not displayed.

## Author

Sebastian Gibb <mail@sebastiangibb.de>

Guillaume Deflandre <guillaume.deflandre@uclouvain.be>

## Examples

``` r

## No modifications
calculateFragments("ACE")
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>          mz ion type pos z seq        peptide
#> 1  72.04439  b1    b   1 1   A AC[+57.02146]E
#> 2 232.07504  b2    b   2 1  AC AC[+57.02146]E
#> 3 148.06043  y1    y   1 1   E AC[+57.02146]E
#> 4 308.09108  y2    y   2 1  CE AC[+57.02146]E
#> 5 130.04987 y1_   y_   1 1   E AC[+57.02146]E
#> 6 290.08052 y2_   y_   2 1  CE AC[+57.02146]E

## Multiple ion types and charge states
calculateFragments("ACE",
                   type = c("a", "b", "c", "x", "y", "z"),
                   z = 1:2)
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>           mz ion type pos z seq        peptide
#> 1   44.04947  a1    a   1 1   A AC[+57.02146]E
#> 2   22.52837  a1    a   1 2   A AC[+57.02146]E
#> 3  204.08012  a2    a   2 1  AC AC[+57.02146]E
#> 4  102.54370  a2    a   2 2  AC AC[+57.02146]E
#> 5   72.04439  b1    b   1 1   A AC[+57.02146]E
#> 6   36.52583  b1    b   1 2   A AC[+57.02146]E
#> 7  232.07504  b2    b   2 1  AC AC[+57.02146]E
#> 8  116.54116  b2    b   2 2  AC AC[+57.02146]E
#> 9   89.07094  c1    c   1 1   A AC[+57.02146]E
#> 10  45.03911  c1    c   1 2   A AC[+57.02146]E
#> 11 249.10158  c2    c   2 1  AC AC[+57.02146]E
#> 12 125.05443  c2    c   2 2  AC AC[+57.02146]E
#> 13 174.03970  x1    x   1 1   E AC[+57.02146]E
#> 14  87.52349  x1    x   1 2   E AC[+57.02146]E
#> 15 334.07035  x2    x   2 1  CE AC[+57.02146]E
#> 16 167.53881  x2    x   2 2  CE AC[+57.02146]E
#> 17 148.06043  y1    y   1 1   E AC[+57.02146]E
#> 18  74.53385  y1    y   1 2   E AC[+57.02146]E
#> 19 308.09108  y2    y   2 1  CE AC[+57.02146]E
#> 20 154.54918  y2    y   2 2  CE AC[+57.02146]E
#> 21 131.03388  z1    z   1 1   E AC[+57.02146]E
#> 22  66.02058  z1    z   1 2   E AC[+57.02146]E
#> 23 291.06453  z2    z   2 1  CE AC[+57.02146]E
#> 24 146.03590  z2    z   2 2  CE AC[+57.02146]E
#> 25 156.02913 x1_   x_   1 1   E AC[+57.02146]E
#> 26  78.51820 x1_   x_   1 2   E AC[+57.02146]E
#> 27 316.05978 x2_   x_   2 1  CE AC[+57.02146]E
#> 28 158.53353 x2_   x_   2 2  CE AC[+57.02146]E
#> 29 130.04987 y1_   y_   1 1   E AC[+57.02146]E
#> 30  65.52857 y1_   y_   1 2   E AC[+57.02146]E
#> 31 290.08052 y2_   y_   2 1  CE AC[+57.02146]E
#> 32 145.54390 y2_   y_   2 2  CE AC[+57.02146]E
#> 33 113.02332 z1_   z_   1 1   E AC[+57.02146]E
#> 34  57.01530 z1_   z_   1 2   E AC[+57.02146]E
#> 35 273.05397 z2_   z_   2 1  CE AC[+57.02146]E
#> 36 137.03062 z2_   z_   2 2  CE AC[+57.02146]E

## Positional modification written directly in the sequence string
## The annotation style must be supported by PTMods::convertAnnotation
calculateFragments("T[+79.966]CE")
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>         mz ion type pos z seq                 peptide
#> 1 182.0210  b1    b   1 1   T T[+79.966]C[+57.02146]E
#> 2 342.0516  b2    b   2 1  TC T[+79.966]C[+57.02146]E
#> 3 148.0604  y1    y   1 1   E T[+79.966]C[+57.02146]E
#> 4 308.0911  y2    y   2 1  CE T[+79.966]C[+57.02146]E
#> 5 130.0499 y1_   y_   1 1   E T[+79.966]C[+57.02146]E
#> 6 290.0805 y2_   y_   2 1  CE T[+79.966]C[+57.02146]E
calculateFragments("T[Phospho]CE")
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>         mz ion type pos z seq                 peptide
#> 1 182.0213  b1    b   1 1   T T[Phospho]C[+57.02146]E
#> 2 342.0519  b2    b   2 1  TC T[Phospho]C[+57.02146]E
#> 3 148.0604  y1    y   1 1   E T[Phospho]C[+57.02146]E
#> 4 308.0911  y2    y   2 1  CE T[Phospho]C[+57.02146]E
#> 5 130.0499 y1_   y_   1 1   E T[Phospho]C[+57.02146]E
#> 6 290.0805 y2_   y_   2 1  CE T[Phospho]C[+57.02146]E
## Notice carbamidomethylation applied by default, but ignored if already
## present.
calculateFragments("T[UNIMOD:21]C[Carbamidomethyl]E")
#>         mz ion type pos z seq                         peptide
#> 1 182.0213  b1    b   1 1   T T[UNIMOD:21]C[Carbamidomethyl]E
#> 2 342.0519  b2    b   2 1  TC T[UNIMOD:21]C[Carbamidomethyl]E
#> 3 148.0604  y1    y   1 1   E T[UNIMOD:21]C[Carbamidomethyl]E
#> 4 308.0911  y2    y   2 1  CE T[UNIMOD:21]C[Carbamidomethyl]E
#> 5 130.0499 y1_   y_   1 1   E T[UNIMOD:21]C[Carbamidomethyl]E
#> 6 290.0805 y2_   y_   2 1  CE T[UNIMOD:21]C[Carbamidomethyl]E

## neutral loss
defaultNeutralLoss()
#> $water
#> [1] "Cterm" "D"     "E"     "S"     "T"    
#> 
#> $ammonia
#> [1] "K" "N" "Q" "R"
#> 

## disable water loss on the C terminal
defaultNeutralLoss(disableWaterLoss="Cterm")
#> $water
#> [1] "D" "E" "S" "T"
#> 
#> $ammonia
#> [1] "K" "N" "Q" "R"
#> 

## real example
calculateFragments("PQR")
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>          mz ion type pos z seq peptide
#> 1  98.06004  b1    b   1 1   P     PQR
#> 2 226.11862  b2    b   2 1  PQ     PQR
#> 3 175.11895  y1    y   1 1   R     PQR
#> 4 303.17753  y2    y   2 1  QR     PQR
#> 5 157.10839 y1_   y_   1 1   R     PQR
#> 6 285.16697 y2_   y_   2 1  QR     PQR
#> 7 286.15098 y2*   y*   2 1  QR     PQR
calculateFragments("PQR",
                   neutralLoss=defaultNeutralLoss(disableWaterLoss="Cterm"))
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>          mz ion type pos z seq peptide
#> 1  98.06004  b1    b   1 1   P     PQR
#> 2 226.11862  b2    b   2 1  PQ     PQR
#> 3 175.11895  y1    y   1 1   R     PQR
#> 4 303.17753  y2    y   2 1  QR     PQR
#> 5 286.15098 y2*   y*   2 1  QR     PQR
calculateFragments("PQR",
                   neutralLoss=defaultNeutralLoss(disableAmmoniaLoss="Q"))
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>          mz ion type pos z seq peptide
#> 1  98.06004  b1    b   1 1   P     PQR
#> 2 226.11862  b2    b   2 1  PQ     PQR
#> 3 175.11895  y1    y   1 1   R     PQR
#> 4 303.17753  y2    y   2 1  QR     PQR
#> 5 157.10839 y1_   y_   1 1   R     PQR
#> 6 285.16697 y2_   y_   2 1  QR     PQR

## disable neutral loss completely
calculateFragments("PQR", neutralLoss=NULL)
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>          mz ion type pos z seq peptide
#> 1  98.06004  b1    b   1 1   P     PQR
#> 2 226.11862  b2    b   2 1  PQ     PQR
#> 3 175.11895  y1    y   1 1   R     PQR
#> 4 303.17753  y2    y   2 1  QR     PQR

## Recommended workflow: use PTMods functions to produce positional sequences
## before calling calculateFragments.

## Fixed modification (Carbamidomethyl on C) using addFixedModifications
seq_fixed <- PTMods::addFixedModifications("ACE",
                                           fixedModifications = c(C = 57.02))
calculateFragments(seq_fixed)
#>          mz ion type pos z seq     peptide
#> 1  72.04439  b1    b   1 1   A AC[+57.02]E
#> 2 232.07358  b2    b   2 1  AC AC[+57.02]E
#> 3 148.06043  y1    y   1 1   E AC[+57.02]E
#> 4 308.08962  y2    y   2 1  CE AC[+57.02]E
#> 5 130.04987 y1_   y_   1 1   E AC[+57.02]E
#> 6 290.07906 y2_   y_   2 1  CE AC[+57.02]E

## Fixed modification including N-terminus using addFixedModifications
seq_nterm <- PTMods::addFixedModifications(
    "ACE",
    fixedModifications = c(C = 57.02, Nterm = 229.16))
calculateFragments(seq_nterm)
#>         mz ion type pos z seq               peptide
#> 1 301.2044  b1    b   1 1   A [+229.16]-AC[+57.02]E
#> 2 461.2336  b2    b   2 1  AC [+229.16]-AC[+57.02]E
#> 3 148.0604  y1    y   1 1   E [+229.16]-AC[+57.02]E
#> 4 308.0896  y2    y   2 1  CE [+229.16]-AC[+57.02]E
#> 5 130.0499 y1_   y_   1 1   E [+229.16]-AC[+57.02]E
#> 6 290.0791 y2_   y_   2 1  CE [+229.16]-AC[+57.02]E

## Variable modification (delta mass on A) using addVariableModifications
seq_var <- PTMods::addVariableModifications("ACE",
                                            variableModifications = c(A = 43.25))
calculateFragments(seq_var)
#> Warning: Could not find Unimod entry for modification +43.25, see `?modifications`
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>           mz ion type pos z seq                peptide
#> 1   72.04439  b1    b   1 1   A         AC[+57.02146]E
#> 2  232.07504  b2    b   2 1  AC         AC[+57.02146]E
#> 3  148.06043  y1    y   1 1   E         AC[+57.02146]E
#> 4  308.09108  y2    y   2 1  CE         AC[+57.02146]E
#> 5  130.04987 y1_   y_   1 1   E         AC[+57.02146]E
#> 6  290.08052 y2_   y_   2 1  CE         AC[+57.02146]E
#> 7  115.29439  b1    b   1 1   A A[+43.25]C[+57.02146]E
#> 8  275.32504  b2    b   2 1  AC A[+43.25]C[+57.02146]E
#> 9  148.06043  y1    y   1 1   E A[+43.25]C[+57.02146]E
#> 10 308.09108  y2    y   2 1  CE A[+43.25]C[+57.02146]E
#> 11 130.04987 y1_   y_   1 1   E A[+43.25]C[+57.02146]E
#> 12 290.08052 y2_   y_   2 1  CE A[+43.25]C[+57.02146]E

## Both fixed and variable modifications using addModifications
seq_mods <- PTMods::addModifications("ARGSHKATC",
                                     fixedModifications = c(C = 57),
                                     variableModifications = c(S = 79, T = 79),
                                     maxMods = 2)
calculateFragments(seq_mods)
#> Warning: Could not find Unimod entry for modification +57, see `?modifications`
#> Warning: Could not find Unimod entry for modification +79, see `?modifications`
#> Warning: Could not find Unimod entry for modification +57, see `?modifications`
#> Warning: Could not find Unimod entry for modification +79, see `?modifications`
#> Warning: Could not find Unimod entry for modification +57, see `?modifications`
#> Fixed modifications used: C=57.02146
#> Variable modifications used: None
#>             mz ion type pos z      seq                        peptide
#> 1     72.04439  b1    b   1 1        A      ARGSHKATC[+57][+57.02146]
#> 2    228.14550  b2    b   2 1       AR      ARGSHKATC[+57][+57.02146]
#> 3    285.16696  b3    b   3 1      ARG      ARGSHKATC[+57][+57.02146]
#> 4    372.19899  b4    b   4 1     ARGS      ARGSHKATC[+57][+57.02146]
#> 5    509.25790  b5    b   5 1    ARGSH      ARGSHKATC[+57][+57.02146]
#> 6    637.35286  b6    b   6 1   ARGSHK      ARGSHKATC[+57][+57.02146]
#> 7    708.38997  b7    b   7 1  ARGSHKA      ARGSHKATC[+57][+57.02146]
#> 8    809.43765  b8    b   8 1 ARGSHKAT      ARGSHKATC[+57][+57.02146]
#> 9    236.04849  y1    y   1 1        C      ARGSHKATC[+57][+57.02146]
#> 10   337.09617  y2    y   2 1       TC      ARGSHKATC[+57][+57.02146]
#> 11   408.13328  y3    y   3 1      ATC      ARGSHKATC[+57][+57.02146]
#> 12   536.22824  y4    y   4 1     KATC      ARGSHKATC[+57][+57.02146]
#> 13   673.28715  y5    y   5 1    HKATC      ARGSHKATC[+57][+57.02146]
#> 14   760.31918  y6    y   6 1   SHKATC      ARGSHKATC[+57][+57.02146]
#> 15   817.34064  y7    y   7 1  GSHKATC      ARGSHKATC[+57][+57.02146]
#> 16   973.44175  y8    y   8 1 RGSHKATC      ARGSHKATC[+57][+57.02146]
#> 17   491.24733 b5_   b_   5 1    ARGSH      ARGSHKATC[+57][+57.02146]
#> 18   619.34229 b6_   b_   6 1   ARGSHK      ARGSHKATC[+57][+57.02146]
#> 19   690.37940 b7_   b_   7 1  ARGSHKA      ARGSHKATC[+57][+57.02146]
#> 20   791.42708 b8_   b_   8 1 ARGSHKAT      ARGSHKATC[+57][+57.02146]
#> 21   390.12272 y3_   y_   3 1      ATC      ARGSHKATC[+57][+57.02146]
#> 22   518.21768 y4_   y_   4 1     KATC      ARGSHKATC[+57][+57.02146]
#> 23   655.27659 y5_   y_   5 1    HKATC      ARGSHKATC[+57][+57.02146]
#> 24   742.30862 y6_   y_   6 1   SHKATC      ARGSHKATC[+57][+57.02146]
#> 25   799.33008 y7_   y_   7 1  GSHKATC      ARGSHKATC[+57][+57.02146]
#> 26   955.43119 y8_   y_   8 1 RGSHKATC      ARGSHKATC[+57][+57.02146]
#> 27   218.03793 y1_   y_   1 1        C      ARGSHKATC[+57][+57.02146]
#> 28   319.08561 y2_   y_   2 1       TC      ARGSHKATC[+57][+57.02146]
#> 29   268.14041 b3*   b*   3 1      ARG      ARGSHKATC[+57][+57.02146]
#> 30   355.17244 b4*   b*   4 1     ARGS      ARGSHKATC[+57][+57.02146]
#> 31   492.23135 b5*   b*   5 1    ARGSH      ARGSHKATC[+57][+57.02146]
#> 32   620.32631 b6*   b*   6 1   ARGSHK      ARGSHKATC[+57][+57.02146]
#> 33   691.36342 b7*   b*   7 1  ARGSHKA      ARGSHKATC[+57][+57.02146]
#> 34   792.41110 b8*   b*   8 1 ARGSHKAT      ARGSHKATC[+57][+57.02146]
#> 35   519.20169 y4*   y*   4 1     KATC      ARGSHKATC[+57][+57.02146]
#> 36   656.26060 y5*   y*   5 1    HKATC      ARGSHKATC[+57][+57.02146]
#> 37   743.29263 y6*   y*   6 1   SHKATC      ARGSHKATC[+57][+57.02146]
#> 38   800.31409 y7*   y*   7 1  GSHKATC      ARGSHKATC[+57][+57.02146]
#> 39   956.41520 y8*   y*   8 1 RGSHKATC      ARGSHKATC[+57][+57.02146]
#> 40    72.04439  b1    b   1 1        A ARGS[+79]HKATC[+57][+57.02146]
#> 41   228.14550  b2    b   2 1       AR ARGS[+79]HKATC[+57][+57.02146]
#> 42   285.16696  b3    b   3 1      ARG ARGS[+79]HKATC[+57][+57.02146]
#> 43   451.19899  b4    b   4 1     ARGS ARGS[+79]HKATC[+57][+57.02146]
#> 44   588.25790  b5    b   5 1    ARGSH ARGS[+79]HKATC[+57][+57.02146]
#> 45   716.35286  b6    b   6 1   ARGSHK ARGS[+79]HKATC[+57][+57.02146]
#> 46   787.38997  b7    b   7 1  ARGSHKA ARGS[+79]HKATC[+57][+57.02146]
#> 47   888.43765  b8    b   8 1 ARGSHKAT ARGS[+79]HKATC[+57][+57.02146]
#> 48   236.04849  y1    y   1 1        C ARGS[+79]HKATC[+57][+57.02146]
#> 49   337.09617  y2    y   2 1       TC ARGS[+79]HKATC[+57][+57.02146]
#> 50   408.13328  y3    y   3 1      ATC ARGS[+79]HKATC[+57][+57.02146]
#> 51   536.22824  y4    y   4 1     KATC ARGS[+79]HKATC[+57][+57.02146]
#> 52   673.28715  y5    y   5 1    HKATC ARGS[+79]HKATC[+57][+57.02146]
#> 53   839.31918  y6    y   6 1   SHKATC ARGS[+79]HKATC[+57][+57.02146]
#> 54   896.34064  y7    y   7 1  GSHKATC ARGS[+79]HKATC[+57][+57.02146]
#> 55  1052.44175  y8    y   8 1 RGSHKATC ARGS[+79]HKATC[+57][+57.02146]
#> 56   570.24733 b5_   b_   5 1    ARGSH ARGS[+79]HKATC[+57][+57.02146]
#> 57   698.34229 b6_   b_   6 1   ARGSHK ARGS[+79]HKATC[+57][+57.02146]
#> 58   769.37940 b7_   b_   7 1  ARGSHKA ARGS[+79]HKATC[+57][+57.02146]
#> 59   870.42708 b8_   b_   8 1 ARGSHKAT ARGS[+79]HKATC[+57][+57.02146]
#> 60   390.12272 y3_   y_   3 1      ATC ARGS[+79]HKATC[+57][+57.02146]
#> 61   518.21768 y4_   y_   4 1     KATC ARGS[+79]HKATC[+57][+57.02146]
#> 62   655.27659 y5_   y_   5 1    HKATC ARGS[+79]HKATC[+57][+57.02146]
#> 63   821.30862 y6_   y_   6 1   SHKATC ARGS[+79]HKATC[+57][+57.02146]
#> 64   878.33008 y7_   y_   7 1  GSHKATC ARGS[+79]HKATC[+57][+57.02146]
#> 65  1034.43119 y8_   y_   8 1 RGSHKATC ARGS[+79]HKATC[+57][+57.02146]
#> 66   218.03793 y1_   y_   1 1        C ARGS[+79]HKATC[+57][+57.02146]
#> 67   319.08561 y2_   y_   2 1       TC ARGS[+79]HKATC[+57][+57.02146]
#> 68   268.14041 b3*   b*   3 1      ARG ARGS[+79]HKATC[+57][+57.02146]
#> 69   434.17244 b4*   b*   4 1     ARGS ARGS[+79]HKATC[+57][+57.02146]
#> 70   571.23135 b5*   b*   5 1    ARGSH ARGS[+79]HKATC[+57][+57.02146]
#> 71   699.32631 b6*   b*   6 1   ARGSHK ARGS[+79]HKATC[+57][+57.02146]
#> 72   770.36342 b7*   b*   7 1  ARGSHKA ARGS[+79]HKATC[+57][+57.02146]
#> 73   871.41110 b8*   b*   8 1 ARGSHKAT ARGS[+79]HKATC[+57][+57.02146]
#> 74   519.20169 y4*   y*   4 1     KATC ARGS[+79]HKATC[+57][+57.02146]
#> 75   656.26060 y5*   y*   5 1    HKATC ARGS[+79]HKATC[+57][+57.02146]
#> 76   822.29263 y6*   y*   6 1   SHKATC ARGS[+79]HKATC[+57][+57.02146]
#> 77   879.31409 y7*   y*   7 1  GSHKATC ARGS[+79]HKATC[+57][+57.02146]
#> 78  1035.41520 y8*   y*   8 1 RGSHKATC ARGS[+79]HKATC[+57][+57.02146]
#> 79    72.04439  b1    b   1 1        A ARGSHKAT[+79]C[+57][+57.02146]
#> 80   228.14550  b2    b   2 1       AR ARGSHKAT[+79]C[+57][+57.02146]
#> 81   285.16696  b3    b   3 1      ARG ARGSHKAT[+79]C[+57][+57.02146]
#> 82   372.19899  b4    b   4 1     ARGS ARGSHKAT[+79]C[+57][+57.02146]
#> 83   509.25790  b5    b   5 1    ARGSH ARGSHKAT[+79]C[+57][+57.02146]
#> 84   637.35286  b6    b   6 1   ARGSHK ARGSHKAT[+79]C[+57][+57.02146]
#> 85   708.38997  b7    b   7 1  ARGSHKA ARGSHKAT[+79]C[+57][+57.02146]
#> 86   888.43765  b8    b   8 1 ARGSHKAT ARGSHKAT[+79]C[+57][+57.02146]
#> 87   236.04849  y1    y   1 1        C ARGSHKAT[+79]C[+57][+57.02146]
#> 88   416.09617  y2    y   2 1       TC ARGSHKAT[+79]C[+57][+57.02146]
#> 89   487.13328  y3    y   3 1      ATC ARGSHKAT[+79]C[+57][+57.02146]
#> 90   615.22824  y4    y   4 1     KATC ARGSHKAT[+79]C[+57][+57.02146]
#> 91   752.28715  y5    y   5 1    HKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 92   839.31918  y6    y   6 1   SHKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 93   896.34064  y7    y   7 1  GSHKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 94  1052.44175  y8    y   8 1 RGSHKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 95   491.24733 b5_   b_   5 1    ARGSH ARGSHKAT[+79]C[+57][+57.02146]
#> 96   619.34229 b6_   b_   6 1   ARGSHK ARGSHKAT[+79]C[+57][+57.02146]
#> 97   690.37940 b7_   b_   7 1  ARGSHKA ARGSHKAT[+79]C[+57][+57.02146]
#> 98   870.42708 b8_   b_   8 1 ARGSHKAT ARGSHKAT[+79]C[+57][+57.02146]
#> 99   469.12272 y3_   y_   3 1      ATC ARGSHKAT[+79]C[+57][+57.02146]
#> 100  597.21768 y4_   y_   4 1     KATC ARGSHKAT[+79]C[+57][+57.02146]
#> 101  734.27659 y5_   y_   5 1    HKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 102  821.30862 y6_   y_   6 1   SHKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 103  878.33008 y7_   y_   7 1  GSHKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 104 1034.43119 y8_   y_   8 1 RGSHKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 105  218.03793 y1_   y_   1 1        C ARGSHKAT[+79]C[+57][+57.02146]
#> 106  398.08561 y2_   y_   2 1       TC ARGSHKAT[+79]C[+57][+57.02146]
#> 107  268.14041 b3*   b*   3 1      ARG ARGSHKAT[+79]C[+57][+57.02146]
#> 108  355.17244 b4*   b*   4 1     ARGS ARGSHKAT[+79]C[+57][+57.02146]
#> 109  492.23135 b5*   b*   5 1    ARGSH ARGSHKAT[+79]C[+57][+57.02146]
#> 110  620.32631 b6*   b*   6 1   ARGSHK ARGSHKAT[+79]C[+57][+57.02146]
#> 111  691.36342 b7*   b*   7 1  ARGSHKA ARGSHKAT[+79]C[+57][+57.02146]
#> 112  871.41110 b8*   b*   8 1 ARGSHKAT ARGSHKAT[+79]C[+57][+57.02146]
#> 113  598.20169 y4*   y*   4 1     KATC ARGSHKAT[+79]C[+57][+57.02146]
#> 114  735.26060 y5*   y*   5 1    HKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 115  822.29263 y6*   y*   6 1   SHKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 116  879.31409 y7*   y*   7 1  GSHKATC ARGSHKAT[+79]C[+57][+57.02146]
#> 117 1035.41520 y8*   y*   8 1 RGSHKATC ARGSHKAT[+79]C[+57][+57.02146]
```

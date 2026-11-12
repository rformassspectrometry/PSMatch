# A class for peptide-spectrum matches

The `PSM` class is a simple class to store and manipulate
peptide-spectrum matches. The class encapsulates PSM data as a
`DataFrame` (or more specifically a `DFrame`) with additional
lightweight metadata annotation.

There are two types of `PSM` objects:

- Objects with duplicated spectrum identifiers. This holds for multiple
  matches to the same spectrum, be it different peptide sequences or the
  same sequence with or without a post-translational modification. Such
  objects are typically created with the `PSM()` constructor starting
  from `mzIdentML` files.

- Reduced objects where the spectrum identifiers (or any equivalent
  column) are unique keys within the PSM table. Matches to the same
  scan/spectrum are merged into a single PSM data row. Reduced `PSM`
  object are created with the `reducePSMs()` function. See examples
  below.

Objects can be checked for their reduced state with the `reduced()`
function which returns `TRUE` for reduced instances, `FALSE` when the
spectrum identifiers are duplicated, or NA when unknown. The flag can
also be set explicitly with the `reduced()<-` setter.

## Usage

``` r
PSM(
  x,
  spectrum = NA,
  peptide = NA,
  protein = NA,
  decoy = NA,
  rank = NA,
  score = NA,
  fdr = NA,
  parser = c("mzR", "mzID"),
  BPPARAM = SerialParam()
)

reduced(object, spectrum = psmVariables(object)["spectrum"])

reduced(object) <- value

psmVariables(object, which = "all")

reducePSMs(object, k = object[[psmVariables(object)["spectrum"]]])

# S4 method for class 'PSM'
adjacencyMatrix(object)
```

## Arguments

- x:

  [`character()`](https://rdrr.io/r/base/character.html) of mzid file
  names, an instance of class `PSM`, or a `data.frame`.

- spectrum:

  `character(1)` variable name that defines a spectrum in the PSM data.
  Default are `"spectrumID"` (mzR parser) or `"spectrumid"` (mzID
  parser). It is also used to calculate the reduced state.

- peptide:

  `character(1)` variable name that defines a peptide in the PSM data.
  Detaults are `"sequence"` (mzR parser) or `"pepSeq"` (mzID parser).

- protein:

  `character(1)` variable name that defines a protein in the PSM data.
  Detaults are `"DatabaseAccess"` (mzR parser) or `"accession"` (mzID
  parser).

- decoy:

  `character(1)` variable name that defines a decoy hit in the PSM data.
  Detaults are `"isDecoy"` (mzR parser) or `"isdecoy"` (mzID parser).

- rank:

  `character(1)` variable name that defines the rank of the peptide
  spectrum match in the PSM data. Default is `"rank"`.

- score:

  `character(1)` variable name that defines the PSM score. This value
  isn't set by default as it depends on the search engine and
  application. Default is `NA`.

- fdr:

  `character(1)` variable name that defines that defines the spectrum
  FDR (or any similar/relevant metric that can be used for filtering).
  This value isn't set by default as it depends on the search engine and
  application. Default is `NA`.

- parser:

  `character(1)` defining the parser to be used to read the `mzIdentML`
  files. One of `"mzR"` (default) or `"mzID"`.

- BPPARAM:

  an object from the `BiocParallel` package to control parallel
  processing. The default value is SerialParam()\` to read files in
  series.

- object:

  An instance of class `PSM`.

- value:

  new value to be passed to setter.

- which:

  [`character()`](https://rdrr.io/r/base/character.html) with the PSM
  variable name to retrieve. If `"all"` (default), all named variables
  are returned. See `PSM()` for valid PSM variables.

- k:

  A `vector` or `factor` of length equal to `nrow(x)` that defines the
  primary key used to reduce `x`. This typically corresponds to the
  spectrum identifier. The defauls is to use the spectrum PSM variable.

## Value

`PSM()` returns a `PSM` object.

`reducePSMs()` returns a reduced version of the `x` input.

## Creating and using PSM objects

- The `PSM()` constructor uses parsers provided by the `mzR` or `mzID`
  packages to read the `mzIdentML` data. The vignette describes some
  apparent differences in their outputs. The constructor input is a
  character of one more multiple file names.

- `PSM` objects can also be created from a `data.frame` object (or any
  variable that can be coerced into a
  [S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html).

- Finally, `PSM()` can also take a `PSM` object as input, which leaves
  the PSM data as is and is used to set/update the PSM variables.

- The constructor can also initialise variables (called *PSM variables*)
  needed for downstream processing, notably filtering (see
  [`filterPSMs()`](https://rformassspectrometry.github.io/PSMatch/reference/filterPSMs.md))
  and to generate a peptide-by-protein adjacency matrix (see
  [`makeAdjacencyMatrix()`](https://rformassspectrometry.github.io/PSMatch/reference/adjacencyMatrix.md)).
  These variables can be extracted with the `psmVariables()` function.
  They represent the columns in the PSM table that identify spectra,
  peptides, proteins, decoy peptides hit ranks and, optionally, a PSM
  score. The value of these variables will depend on the backend used to
  create the object, or left blank (i.e. encoded as `NA`) when building
  an object by hand from a `data.frame`. In such situation, they need to
  be passed explicitly by the user as arguments to `PSM()`.

- The
  [`adjacencyMatrix()`](https://rformassspectrometry.github.io/PSMatch/reference/adjacencyMatrix.md)
  accessor can be used to retrieve the binary sparse peptide-by-protein
  adjacency matrix from the PSM object. It also relies on PSM variables
  which thus need to be set beforehand. For more flexibility in the
  generation of the adjacency matrix (for non-binary matrices), use
  [`makeAdjacencyMatrix()`](https://rformassspectrometry.github.io/PSMatch/reference/adjacencyMatrix.md).

## Examples

``` r

## ---------------------------------
## Example with a single mzid file
## ---------------------------------

f <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache

## mzR parser (default)
psm <- PSM(f)
psm
#> PSM with 5802 rows and 35 columns.
#> names(35): sequence spectrumID ... subReplacementResidue subLocation

## PSM variables
psmVariables(psm)
#>         spectrum          peptide          protein            decoy 
#>     "spectrumID"       "sequence" "DatabaseAccess"        "isDecoy" 
#>             rank            score              fdr 
#>           "rank"               NA               NA 

## mzID parser
psm_mzid <- PSM(f, parser = "mzID")
#> Loading required namespace: mzID
#> reading 5b71a068abd_7857...
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#> Warning: 'as.is' should be specified by the caller; using TRUE
#>  DONE!
psm_mzid
#> PSM with 5759 rows and 30 columns.
#> names(30): spectrumid scan number(s) ... spectrumFile databaseFile

## different PSM variables
psmVariables(psm_mzid)
#>     spectrum      peptide      protein        decoy         rank        score 
#> "spectrumid"     "pepseq"  "accession"    "isdecoy"       "rank"           NA 
#>          fdr 
#>           NA 

## Reducing the PSM data
(i <- which(duplicated(psm$spectrumID))[1:2])
#> [1]  9 24
(i <- which(psm$spectrumID %in% psm$spectrumID[i]))
#> [1]  8  9 23 24 25 26
psm2 <- psm[i, ]
reduced(psm2)
#> [1] NA

## Peptide sequence CIDRARHVEVQIFGDGKGRVVALGERDCSLQRR with
## Carbamidomethyl modifications at positions 1 and 28.
DataFrame(psm2[, c("sequence", "spectrumID", "modName", "modLocation")])
#> DataFrame with 6 rows and 4 columns
#>                sequence    spectrumID       modName modLocation
#>             <character>   <character>   <character>   <integer>
#> EH7807.8  CIDRARHVEV... controller... Carbamidom...           1
#> EH7807.9  CIDRARHVEV... controller... Carbamidom...          28
#> EH7807.23 KCNQCLKVAC... controller... Carbamidom...           2
#> EH7807.24 KCNQCLKVAC... controller... Carbamidom...           5
#> EH7807.25 KCNQCLKVAC... controller... Carbamidom...          10
#> EH7807.26 KCNQCLKVAC... controller... Carbamidom...          15
reduced(psm2) <- FALSE
reduced(psm2)
#> [1] FALSE

## uses by default the spectrum PSM variable, as defined during
## the construction - see psmVariables()
rpsm2 <- reducePSMs(psm2)
rpsm2
#> Reduced PSM with 2 rows and 35 columns.
#> names(35): sequence spectrumID ... subReplacementResidue subLocation
DataFrame(rpsm2[, c("sequence", "spectrumID", "modName", "modLocation")])
#> DataFrame with 2 rows and 4 columns
#>                                                    sequence    spectrumID
#>                                                 <character>   <character>
#> controllerType=0 controllerNumber=1 scan=5490 KCNQCLKVAC... controller...
#> controllerType=0 controllerNumber=1 scan=6518 CIDRARHVEV... controller...
#>                                                     modName   modLocation
#>                                                 <character> <IntegerList>
#> controllerType=0 controllerNumber=1 scan=5490 Carbamidom...    2,5,10,...
#> controllerType=0 controllerNumber=1 scan=6518 Carbamidom...          1,28
reduced(rpsm2)
#> [1] TRUE

## ---------------------------------
## Multiple mzid files
## ---------------------------------

library(rpx)
PXD022816 <- PXDataset("PXD022816")
#> Loading PXD022816 from cache.
PXD022816
#> Project PXD022816 with 33 files
#>  
#> Resource ID BFC1 in cache in /github/home/.cache/R/rpx.
#>  [1] 'PXD022816_community_annotated.sdrf.tsv' ... [33] 'checksum.txt'
#>  Use 'pxfiles(.)' to see all files.

(mzids <- pxget(PXD022816, grep("mzID", pxfiles(PXD022816))[1:2]))
#> Project PXD022816 files (33):
#>  [remote] PXD022816_community_annotated.sdrf.tsv
#>  [local]  QEP2LC6_HeLa_50ng_251120_01-calib.mzID.gz
#>  [remote] QEP2LC6_HeLa_50ng_251120_01-calib.mzML
#>  [remote] QEP2LC6_HeLa_50ng_251120_01.raw
#>  [local]  QEP2LC6_HeLa_50ng_251120_02-calib.mzID.gz
#>  [remote] QEP2LC6_HeLa_50ng_251120_02-calib.mzML
#>  [remote] QEP2LC6_HeLa_50ng_251120_02.raw
#>  [remote] QEP2LC6_HeLa_50ng_251120_03-calib.mzID.gz
#>  [remote] QEP2LC6_HeLa_50ng_251120_03-calib.mzML
#>  [remote] QEP2LC6_HeLa_50ng_251120_03.raw
#>  ...
#> Loading QEP2LC6_HeLa_50ng_251120_01-calib.mzID.gz from cache.
#> Loading QEP2LC6_HeLa_50ng_251120_02-calib.mzID.gz from cache.
#> [1] "/github/home/.cache/R/rpx/f9e5b0a3c61_QEP2LC6_HeLa_50ng_251120_01-calib.mzID.gz"
#> [2] "/github/home/.cache/R/rpx/f9e347dab22_QEP2LC6_HeLa_50ng_251120_02-calib.mzID.gz"
psm <- PSM(mzids)
psm
#> PSM with 50314 rows and 31 columns.
#> names(31): sequence spectrumID ... subReplacementResidue subLocation
psmVariables(psm)
#>         spectrum          peptide          protein            decoy 
#>     "spectrumID"       "sequence" "DatabaseAccess"        "isDecoy" 
#>             rank            score              fdr 
#>           "rank"               NA               NA 

## Here, spectrum identifiers are repeated accross files
psm[grep("scan=20000", psm$spectrumID), "spectrumFile"]
#> [1] "P:\\Projects\\Morgenstern\\HF2_QC_for_the_paper\\QEP2\\2020-11-29-13-56-12\\Task1-CalibrateTask\\QEP2LC6_HeLa_50ng_251120_01-calib.mzML"
#> [2] "P:\\Projects\\Morgenstern\\HF2_QC_for_the_paper\\QEP2\\2020-11-29-13-56-12\\Task1-CalibrateTask\\QEP2LC6_HeLa_50ng_251120_02-calib.mzML"

## Let's create a new primary identifier composed of the scan
## number and the file name
psm$pkey <- paste(sub("^.+Task\\\\", "", psm$spectrumFile),
                  sub("^.+scan=", "", psm$spectrumID),
                  sep = "::")
head(psm$pkey)
#> [1] "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::3426" 
#> [2] "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::20165"
#> [3] "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::20180"
#> [4] "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::15180"
#> [5] "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::10327"
#> [6] "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::12894"

## the PSM is not reduced
reduced(psm, "pkey")
#> [1] FALSE
DataFrame(psm[6:7, ])
#> DataFrame with 2 rows and 32 columns
#>        sequence    spectrumID chargeState      rank passThreshold
#>     <character>   <character>   <integer> <integer>     <logical>
#> 1 SEGSSCALES... controller...           4         1          TRUE
#> 2 SEGSSCALES... controller...           4         1          TRUE
#>   experimentalMassToCharge calculatedMassToCharge  peptideRef    modNum
#>                  <numeric>              <numeric> <character> <integer>
#> 1                  796.367                796.367        P_19         2
#> 2                  796.367                796.367        P_19         2
#>     isDecoy        post         pre     start       end DatabaseAccess
#>   <logical> <character> <character> <integer> <integer>    <character>
#> 1     FALSE           I           K       220       251         Q58D62
#> 2     FALSE           I           K       220       251         Q58D62
#>   DBseqLength   DatabaseSeq DatabaseDescription scan.start.time acquisitionNum
#>     <integer>   <character>         <character>     <character>      <numeric>
#> 1         387 MNVLLLLVLC...       Q58D62|FET...       32.788864          12894
#> 2         387 MNVLLLLVLC...       Q58D62|FET...       32.788864          12894
#>     spectrumFile        idFile MetaMorpheus.score PSM.level.q.value
#>      <character>   <character>          <numeric>         <numeric>
#> 1 P:\\Project... f9e5b0a3c6...            31.2586                 0
#> 2 P:\\Project... f9e5b0a3c6...            31.2586                 0
#>   modPeptideRef       modName   modMass modLocation subOriginalResidue
#>     <character>   <character> <numeric>   <integer>        <character>
#> 1          P_19 Carbamidom...   57.0215           6                 NA
#> 2          P_19 Carbamidom...   57.0215          19                 NA
#>   subReplacementResidue subLocation          pkey
#>             <character>   <integer>   <character>
#> 1                    NA          NA QEP2LC6_He...
#> 2                    NA          NA QEP2LC6_He...

## same sequence, same spectrumID, same file
psm$sequence[6:7]
#> [1] "SEGSSCALESPGSVPVGICHGSLGEPQGNQGK" "SEGSSCALESPGSVPVGICHGSLGEPQGNQGK"
psm$pkey[6:7]
#> [1] "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::12894"
#> [2] "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::12894"

## different modification locations
psm$modLocation[6:7]
#> [1]  6 19

## here, we need to *explicitly* set pkey to reduce
rpsm <- reducePSMs(psm, psm$pkey)
rpsm
#> Reduced PSM with 32739 rows and 32 columns.
#> names(32): sequence spectrumID ... subLocation pkey
reduced(rpsm, "pkey")
#> [1] TRUE

## the two rows are now merged into a single one; the distinct
## modification locations are preserved.
(i <- which(rpsm$pkey == "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::12894"))
#> [1] 2343
DataFrame(rpsm[i, c("sequence", "pkey", "modName", "modLocation")])
#> DataFrame with 1 row and 4 columns
#>                                                      sequence          pkey
#>                                               <CharacterList>   <character>
#> QEP2LC6_HeLa_50ng_251120_01-calib.mzML::12894   SEGSSCALES... QEP2LC6_He...
#>                                                       modName   modLocation
#>                                               <CharacterList> <IntegerList>
#> QEP2LC6_HeLa_50ng_251120_01-calib.mzML::12894   Carbamidom...          6,19

## ---------------------------------
## PSM from a data.frame
## ---------------------------------

psmdf <- data.frame(spectrum = paste0("sp", 1:10),
                    sequence = replicate(10,
                                         paste(sample(getAminoAcids()[-1, "AA"], 10),
                                               collapse = "")),
                    protein = sample(paste0("Prot", LETTERS[1:7]), 10,
                                     replace = TRUE),
                    decoy = rep(FALSE, 10),
                    rank = rep(1, 10),
                    score = runif(10))
psmdf
#>    spectrum   sequence protein decoy rank      score
#> 1       sp1 NHTPAMIQWR   ProtC FALSE    1 0.95426923
#> 2       sp2 NILRMESGYQ   ProtF FALSE    1 0.59342788
#> 3       sp3 IRDPGQSTAH   ProtF FALSE    1 0.71910708
#> 4       sp4 FVMERSYWNL   ProtA FALSE    1 0.04937636
#> 5       sp5 ANEIGSTYOL   ProtG FALSE    1 0.44901218
#> 6       sp6 NKUHARMVGP   ProtC FALSE    1 0.80746808
#> 7       sp7 CMUWPSGDLK   ProtD FALSE    1 0.64619430
#> 8       sp8 IYATDNSGKW   ProtA FALSE    1 0.47634646
#> 9       sp9 ACLUNIHWKG   ProtC FALSE    1 0.15341002
#> 10     sp10 YWKQUECRSD   ProtA FALSE    1 0.59587423

psm <- PSM(psmdf)
psm
#> PSM with 10 rows and 6 columns.
#> names(6): spectrum sequence ... rank score
psmVariables(psm)
#> spectrum  peptide  protein    decoy     rank    score      fdr 
#>       NA       NA       NA       NA       NA       NA       NA 

## no PSM variables set
try(adjacencyMatrix(psm))
#> Error in .local(object, ...) : 
#>   Please define the 'protein' and 'peptide' PSM variables.

## set PSM variables
psm <- PSM(psm, spectrum = "spectrum", peptide = "sequence",
           protein = "protein", decoy = "decoy", rank = "rank")
psm
#> PSM with 10 rows and 6 columns.
#> names(6): spectrum sequence ... rank score
psmVariables(psm)
#>   spectrum    peptide    protein      decoy       rank      score        fdr 
#> "spectrum" "sequence"  "protein"    "decoy"     "rank"         NA         NA 

adjacencyMatrix(psm)
#> 10 x 5 sparse Matrix of class "dgCMatrix"
#>            ProtC ProtF ProtA ProtG ProtD
#> NHTPAMIQWR     1     .     .     .     .
#> NILRMESGYQ     .     1     .     .     .
#> IRDPGQSTAH     .     1     .     .     .
#> FVMERSYWNL     .     .     1     .     .
#> ANEIGSTYOL     .     .     .     1     .
#> NKUHARMVGP     1     .     .     .     .
#> CMUWPSGDLK     .     .     .     .     1
#> IYATDNSGKW     .     .     1     .     .
#> ACLUNIHWKG     1     .     .     .     .
#> YWKQUECRSD     .     .     1     .     .
```

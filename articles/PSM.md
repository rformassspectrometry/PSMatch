# Working with PSM data

**Package**:
*[PSMatch](https://bioconductor.org/packages/3.23/PSMatch)*\
**Authors**: Laurent Gatto \[aut, cre\] (ORCID:
<https://orcid.org/0000-0002-1520-2268>), Johannes Rainer \[aut\]
(ORCID: <https://orcid.org/0000-0002-6977-7147>), Sebastian Gibb \[aut\]
(ORCID: <https://orcid.org/0000-0001-7406-4443>), Samuel Wieczorek
\[ctb\], Thomas Burger \[ctb\], Guillaume Deflandre \[aut\] (ORCID:
<https://orcid.org/0009-0008-1257-2416>)\
**Last modified:** 2026-11-12 13:10:28.559227\
**Compiled**: Thu Nov 12 13:25:02 2026

## Installation instructions

To install the package from Bioconductor, make sure you have the
`BiocManager` package, available from CRAN, and then run

``` r

BiocManager::install("PSMatch")
```

## Introduction

This vignette is one among several illustrating how to use the `PSMatch`
package, focusing on the handling and processing of proteomics
identification data using the `PSM` class. For a general overview of the
package, see the `PSMatch` package manual page
([`?PSMatch`](https://rformassspectrometry.github.io/PSMatch/reference/PSMatch.md))
and references therein.

## Handling and processing identification data

### Loading PSM data

We are going to use an `mzid` file from the `MsDataHub` package.

``` r

f <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
```

    ## see ?MsDataHub and browseVignettes('MsDataHub') for documentation

    ## loading from cache

The
[`PSM()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md)
function parses one or multiple `mzid` files and returns an object of
class `PSM`. This class is a simple extension of the `DFrame` class,
representing the peptide-spectrum matches in a tabular format.

``` r

library("PSMatch")
id <- PSM(f)
```

    ## Loading required namespace: mzR

``` r

id
```

    ## PSM with 5802 rows and 35 columns.
    ## names(35): sequence spectrumID ... subReplacementResidue subLocation

This table contains 5802 matches for 5343 scans and 4938 peptides
sequences, each annotated by 35 variables.

``` r

nrow(id) ## number of matches
```

    ## [1] 5802

``` r

length(unique(id$spectrumID)) ## number of scans
```

    ## [1] 5343

``` r

length(unique(id$sequence))   ## number of peptide sequences
```

    ## [1] 4938

``` r

names(id)
```

    ##  [1] "sequence"                 "spectrumID"              
    ##  [3] "chargeState"              "rank"                    
    ##  [5] "passThreshold"            "experimentalMassToCharge"
    ##  [7] "calculatedMassToCharge"   "peptideRef"              
    ##  [9] "modNum"                   "isDecoy"                 
    ## [11] "post"                     "pre"                     
    ## [13] "start"                    "end"                     
    ## [15] "DatabaseAccess"           "DBseqLength"             
    ## [17] "DatabaseSeq"              "DatabaseDescription"     
    ## [19] "scan.number.s."           "acquisitionNum"          
    ## [21] "spectrumFile"             "idFile"                  
    ## [23] "MS.GF.RawScore"           "MS.GF.DeNovoScore"       
    ## [25] "MS.GF.SpecEValue"         "MS.GF.EValue"            
    ## [27] "MS.GF.QValue"             "MS.GF.PepQValue"         
    ## [29] "modPeptideRef"            "modName"                 
    ## [31] "modMass"                  "modLocation"             
    ## [33] "subOriginalResidue"       "subReplacementResidue"   
    ## [35] "subLocation"

The PSM data are read as is, without any filtering. As we can see below,
we still have all the hits from the forward and reverse (decoy)
databases.

``` r

table(id$isDecoy)
```

    ## 
    ## FALSE  TRUE 
    ##  2906  2896

### Keeping all matches

The data also contains multiple matches for several spectra. The table
below shows the number of individual MS scans that have 1, 2, … up to 5
matches.

``` r

table(table(id$spectrumID))
```

    ## 
    ##    1    2    3    4    5 
    ## 4936  369   26   10    2

More specifically, we can see below how scan 1774 has 4 matches, all to
sequence `RTRYQAEVR`, which itself matches to 4 different proteins:

``` r

i <- grep("scan=1774", id$spectrumID)
id[i, ]
```

    ## PSM with 4 rows and 35 columns.
    ## names(35): sequence spectrumID ... subReplacementResidue subLocation

``` r

id[i, "DatabaseAccess"]
```

    ## [1] "ECA2104" "ECA2867" "ECA3427" "ECA4142"

If the goal is to keep all the matches, but arranged by scan/spectrum,
one can *reduce* the `DataFrame` object by the `spectrumID` variable, so
that each scan correponds to a single row that still stores all
values[^1]:

``` r

id2 <- reducePSMs(id, id$spectrumID)
rownames(id2) <- NULL ## rownames not needed here
dim(id2)
```

    ## [1] 5343   35

The resulting object contains a single entrie for scan 1774 with
information for the multiple matches stored as a list within the table
cell.

``` r

j <- grep("scan=1774", id2$spectrumID)
id2[j, ]
```

    ## Reduced PSM with 1 rows and 35 columns.
    ## names(35): sequence spectrumID ... subReplacementResidue subLocation

``` r

id2[j, "DatabaseAccess"]
```

    ## CharacterList of length 1
    ## [["controllerType=0 controllerNumber=1 scan=1774"]] ECA2104 ECA2867 ECA3427 ECA4142

The identification data could be used to annotate an raw mass
spectrometry `Spectra` object (see the
[`Spectra::joinSpectraData()`](https://rdrr.io/pkg/Spectra/man/combineSpectra.html)
function for details).

### Filtering data

Often, the PSM data is filtered to only retain reliable matches. The
`MSnID` package can be used to set thresholds to attain user-defined
PSM, peptide or protein-level FDRs. Here, we will simply filter out
wrong or the least reliable identifications.

#### Remove decoy hits

``` r

id <- filterPsmDecoy(id)
```

    ## Removed 2896 decoy hits.

``` r

id
```

    ## PSM with 2906 rows and 35 columns.
    ## names(35): sequence spectrumID ... subReplacementResidue subLocation

#### Keep first rank matches

``` r

id <- filterPsmRank(id)
```

    ## Removed 155 PSMs with rank > 1.

``` r

id
```

    ## PSM with 2751 rows and 35 columns.
    ## names(35): sequence spectrumID ... subReplacementResidue subLocation

#### Remove shared peptides

The data still contains shared peptides, i.e. those that different
proteins. For example `QKTRCATRAFKANKGRAR` matches proteins `ECA2869`
and `ECA4278`.

``` r

id[id$sequence == "QKTRCATRAFKANKGRAR", "DatabaseAccess"]
```

    ## [1] "ECA2869" "ECA4278"

We can filter these out to retain unique peptides.

``` r

id <- filterPsmShared(id)
```

    ## Removed 85 shared peptides.

``` r

id
```

    ## PSM with 2666 rows and 35 columns.
    ## names(35): sequence spectrumID ... subReplacementResidue subLocation

This last filtering leaves us with 2666 PSMs.

Note that the `ConnectedComponents` approach defined in this package
allows one to explore protein groups defined by such shared peptides
more thoroughly and make informed decision as to which shared peptides
to retain and which ones to drop. For details see the related vignette
or the
[`ConnectedComponents`](https://rformassspectrometry.github.io/PSMatch/reference/ConnectedComponents.html)
manual page.

#### All filters in one function

This can also be achieved with the
[`filterPSMs()`](https://rformassspectrometry.github.io/PSMatch/reference/filterPSMs.md)
function:

``` r

id <- PSM(f)
filterPSMs(id)
```

    ## Starting with 5802 PSMs:

    ## Removed 2896 decoy hits.

    ## Removed 155 PSMs with rank > 1.

    ## Removed 85 shared peptides.

    ## 2666 PSMs left.

    ## PSM with 2666 rows and 35 columns.
    ## names(35): sequence spectrumID ... subReplacementResidue subLocation

## The `mzR` and `mzID` parsers

The
[`PSM()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md)
function can take two different values for the `parser` parameter,
namely `"mzR"` (the default value) and `"mzID"`.

- **mzR** uses the `openIDfile()` function from the
  *[mzR](https://bioconductor.org/packages/3.23/mzR)* to parse the
  `mzId` file(s), and then coerces the data to a `data.frame` which is
  eventually returned as a `PSM` object. The parser function uses
  dedicated code from the Proteowizard project (included in `mzR`) and
  is generally the fastest approach.

- **mzID** parses the `mzId` file with `mzID()` function from the
  *[mzID](https://bioconductor.org/packages/3.23/mzID)* package, and
  then flattens the data to a `data.frame` with
  [`mzID::flatten()`](https://rdrr.io/pkg/mzID/man/flatten.html) and
  eventuelly returns a `PSM` object. The `mzID` package relies on the
  *[XML](https://CRAN.R-project.org/package=XML)* package. Is is slower
  but is is more robust to variations in the `mzID` implementation, as
  is thus a useful backup when the `mzR` backend fails.

``` r

system.time(id1 <- PSM(f, parser = "mzR"))
```

    ##    user  system elapsed 
    ##   0.167   0.007   0.176

``` r

system.time(id2 <- PSM(f, parser = "mzID"))
```

    ## Loading required namespace: mzID

    ## reading 52023f38f18_7857...

    ##  DONE!

    ##    user  system elapsed 
    ##   5.973   0.066   6.065

Other differences in the two parsers include the columns that are
returned, the way they name them, and, as will shown below the matches
that are returned. Note for instance (and this will be important later),
that there is no equivalent of `"modLocation"` in `id2`.

``` r

names(id1)
```

    ##  [1] "sequence"                 "spectrumID"              
    ##  [3] "chargeState"              "rank"                    
    ##  [5] "passThreshold"            "experimentalMassToCharge"
    ##  [7] "calculatedMassToCharge"   "peptideRef"              
    ##  [9] "modNum"                   "isDecoy"                 
    ## [11] "post"                     "pre"                     
    ## [13] "start"                    "end"                     
    ## [15] "DatabaseAccess"           "DBseqLength"             
    ## [17] "DatabaseSeq"              "DatabaseDescription"     
    ## [19] "scan.number.s."           "acquisitionNum"          
    ## [21] "spectrumFile"             "idFile"                  
    ## [23] "MS.GF.RawScore"           "MS.GF.DeNovoScore"       
    ## [25] "MS.GF.SpecEValue"         "MS.GF.EValue"            
    ## [27] "MS.GF.QValue"             "MS.GF.PepQValue"         
    ## [29] "modPeptideRef"            "modName"                 
    ## [31] "modMass"                  "modLocation"             
    ## [33] "subOriginalResidue"       "subReplacementResidue"   
    ## [35] "subLocation"

``` r

names(id2)
```

    ##  [1] "spectrumid"                "scan number(s)"           
    ##  [3] "acquisitionnum"            "passthreshold"            
    ##  [5] "rank"                      "calculatedmasstocharge"   
    ##  [7] "experimentalmasstocharge"  "chargestate"              
    ##  [9] "ms-gf:denovoscore"         "ms-gf:evalue"             
    ## [11] "ms-gf:pepqvalue"           "ms-gf:qvalue"             
    ## [13] "ms-gf:rawscore"            "ms-gf:specevalue"         
    ## [15] "assumeddissociationmethod" "isotopeerror"             
    ## [17] "isdecoy"                   "post"                     
    ## [19] "pre"                       "end"                      
    ## [21] "start"                     "accession"                
    ## [23] "length"                    "description"              
    ## [25] "pepseq"                    "modified"                 
    ## [27] "modification"              "idFile"                   
    ## [29] "spectrumFile"              "databaseFile"

We also have different number of matches in the two tables:

``` r

nrow(id1)
```

    ## [1] 5802

``` r

nrow(id2)
```

    ## [1] 5759

``` r

table(id1$isDecoy)
```

    ## 
    ## FALSE  TRUE 
    ##  2906  2896

``` r

table(id2$isdecoy)
```

    ## 
    ## FALSE  TRUE 
    ##  2886  2873

Let’s first filter the PSM tables to facilitate focus the comparison of
relevant scans. Note that the default
[`filterPSMs()`](https://rformassspectrometry.github.io/PSMatch/reference/filterPSMs.md)
arguments are set to work with both parser.

``` r

id1_filtered <- filterPSMs(id1)
```

    ## Starting with 5802 PSMs:

    ## Removed 2896 decoy hits.

    ## Removed 155 PSMs with rank > 1.

    ## Removed 85 shared peptides.

    ## 2666 PSMs left.

``` r

id2_filtered <- filterPSMs(id2)
```

    ## Starting with 5759 PSMs:

    ## Removed 2873 decoy hits.

    ## Removed 155 PSMs with rank > 1.

    ## Removed 85 shared peptides.

    ## 2646 PSMs left.

As can be seen, we are also left with 2666 vs 2646 PSMs after filtering.

The difference doesn’t stem from different scans, given that the spectum
identifiers are identical in both tables:

``` r

identical(sort(unique(id1_filtered$spectrumID)),
          sort(unique(id2_filtered$spectrumid)))
```

    ## [1] TRUE

The difference is obvious when we tally a table of spectrum id
occurences in the filtered tables. In `id2_filtered`, each scan is
unique, i.e matched only once.

``` r

anyDuplicated(id2_filtered$spectrumid)
```

    ## [1] 0

However, for `id1_filtered`, we see that some scans are still repeat up
to 4 times in the table:

``` r

table(table(id1_filtered$spectrumID))
```

    ## 
    ##    1    2    3    4 
    ## 2630   13    2    1

The example below shows that these differences stem from the
modification location (`"modLocation"`), that is not report by the
`mzID` parser:

``` r

k <- names(which(table(id1_filtered$spectrumID) == 4))
id1_filtered[id1_filtered$spectrumID == k, "sequence"]
```

    ## [1] "KCNQCLKVACTLFYCK" "KCNQCLKVACTLFYCK" "KCNQCLKVACTLFYCK" "KCNQCLKVACTLFYCK"

``` r

id1_filtered[id1_filtered$spectrumID == k, "modLocation"]
```

    ## [1]  2  5 10 15

``` r

id1_filtered[id1_filtered$spectrumID == k, "modName"]
```

    ## [1] "Carbamidomethyl" "Carbamidomethyl" "Carbamidomethyl" "Carbamidomethyl"

If we remove the `"modLocation"` column, we recoved the same number of
PSMs than with the `mzID` parser.

``` r

id1_filtered$modLocation <- NULL
nrow(unique(id1_filtered))
```

    ## [1] 2646

``` r

nrow(unique(id2_filtered))
```

    ## [1] 2646

## License

This vignette is distributed under a [CC BY-SA
license](https://creativecommons.org/licenses/by-sa/4.0/).

## Session information

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] MsDataHub_1.13.3     PSMatch_1.17.4       PTMods_1.1.1        
    ## [4] S4Vectors_0.51.10    BiocGenerics_0.59.12 generics_0.1.4      
    ## [7] BiocStyle_2.41.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] DBI_1.3.0                   httr2_1.3.0                
    ##   [3] rlang_1.3.0                 magrittr_2.0.5             
    ##   [5] clue_0.3-68                 otel_0.3.0                 
    ##   [7] matrixStats_1.5.0           compiler_4.6.1             
    ##   [9] RSQLite_3.53.3              png_0.1-9                  
    ##  [11] systemfonts_1.3.2           vctrs_0.7.3                
    ##  [13] reshape2_1.4.5              stringr_1.6.0              
    ##  [15] ProtGenerics_1.45.0         pkgconfig_2.0.3            
    ##  [17] MetaboCoreUtils_1.21.2      crayon_1.5.3               
    ##  [19] fastmap_1.2.0               dbplyr_2.6.0               
    ##  [21] XVector_0.53.0              rmarkdown_2.32             
    ##  [23] ragg_1.5.2                  purrr_1.2.2                
    ##  [25] bit_4.6.0                   xfun_0.62                  
    ##  [27] MultiAssayExperiment_1.39.1 cachem_1.1.0               
    ##  [29] jsonlite_2.0.0              blob_1.3.0                 
    ##  [31] DelayedArray_0.39.8         BiocParallel_1.47.0        
    ##  [33] parallel_4.6.1              cluster_2.1.8.3            
    ##  [35] R6_2.6.1                    bslib_0.12.0               
    ##  [37] stringi_1.8.9               GenomicRanges_1.65.4       
    ##  [39] jquerylib_0.1.4             iterators_1.0.14           
    ##  [41] Rcpp_1.1.2                  Seqinfo_1.3.2              
    ##  [43] bookdown_0.48               SummarizedExperiment_1.43.0
    ##  [45] knitr_1.52                  IRanges_2.47.5             
    ##  [47] BiocBaseUtils_1.15.1        Matrix_1.7-6               
    ##  [49] igraph_2.3.4                tidyselect_1.2.1           
    ##  [51] abind_1.4-8                 yaml_2.3.12                
    ##  [53] doParallel_1.0.17           codetools_0.2-20           
    ##  [55] curl_8.1.0                  lattice_0.23-1             
    ##  [57] tibble_3.3.1                plyr_1.8.9                 
    ##  [59] withr_3.0.3                 Biobase_2.73.2             
    ##  [61] KEGGREST_1.53.6             evaluate_1.0.5             
    ##  [63] desc_1.4.3                  Spectra_1.23.7             
    ##  [65] BiocFileCache_3.3.2         ExperimentHub_3.3.2        
    ##  [67] Biostrings_2.81.9           pillar_1.11.1              
    ##  [69] BiocManager_1.30.27         filelock_1.0.3             
    ##  [71] MatrixGenerics_1.25.0       foreach_1.5.2              
    ##  [73] ncdf4_1.24                  BiocVersion_3.23.1         
    ##  [75] glue_1.8.1                  lazyeval_0.2.3             
    ##  [77] tools_4.6.1                 AnnotationHub_4.3.3        
    ##  [79] mzID_1.51.0                 data.table_1.18.6.1        
    ##  [81] mzR_2.47.1                  QFeatures_1.23.2           
    ##  [83] XML_3.99-0.25               fs_2.1.0                   
    ##  [85] grid_4.6.1                  tidyr_1.3.2                
    ##  [87] MsCoreUtils_1.25.4          AnnotationDbi_1.75.2       
    ##  [89] cli_3.6.6                   rappdirs_0.3.4             
    ##  [91] textshaping_1.0.5           S4Arrays_1.13.2            
    ##  [93] dplyr_1.2.1                 AnnotationFilter_1.37.0    
    ##  [95] sass_0.4.10                 digest_0.6.39              
    ##  [97] SparseArray_1.13.4          htmlwidgets_1.6.4          
    ##  [99] memoise_2.0.1               htmltools_0.5.9            
    ## [101] pkgdown_2.2.1.9000          lifecycle_1.0.5            
    ## [103] httr_1.4.9                  bit64_4.8.6                
    ## [105] MASS_7.3-66

[^1]: The rownames aren’t needed here and are removed to reduce to
    output in the the next code chunks displaying parts of `id2`.

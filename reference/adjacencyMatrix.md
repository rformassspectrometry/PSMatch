# Convert to/from an adjacency matrix.

There are two ways that peptide/protein matches are commonly stored:
either as a vector or an adjacency matrix. The functions described below
convert between these two format.

## Usage

``` r
makeAdjacencyMatrix(
  x,
  split = ";",
  peptide = psmVariables(x)["peptide"],
  protein = psmVariables(x)["protein"],
  score = psmVariables(x)["score"],
  binary = FALSE
)

makePeptideProteinVector(m, collapse = ";")

plotAdjacencyMatrix(
  m,
  protColors = 0,
  pepColors = NULL,
  layout = igraph::layout_nicely
)
```

## Arguments

- x:

  Either an instance of class `PSM` or a `character`. See example below
  for details.

- split:

  `character(1)` defining how to split the string of protein identifiers
  (using [`strsplit()`](https://rdrr.io/r/base/strsplit.html)). Default
  is ";". If `NULL`, splitting is ignored.

- peptide:

  `character(1)` indicating the name of the variable that defines
  peptides in the `PSM` object. Default is the `peptide` PSM variable as
  defined in
  [`psmVariables()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).

- protein:

  `character(1)` indicating the name of the variable that defines
  proteins in the `PSM` object. Default is the `peptide` PSM variable as
  defined in
  [`psmVariables()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).

- score:

  `character(1)` indicating the name of the variable that defines PSM
  scores in the `PSM` object. Default is the `score` PSM variable as
  defined in
  [`psmVariables()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).
  Ignored when `NA` (which is the default value unless set by the user
  when constructing the `PSM` object).

- binary:

  `logical(1)` indicates if the adjacency matrix should be strictly
  binary. In such a case, PSMs matching the same peptide but from
  different precursors (for example charge 2 and 3) or carrying
  different PTMs, are counted only once. Default if `FALSE`. This also
  overrides any `score` that would be set.

- m:

  A peptide-by-protein adjacency matrix.

- collapse:

  `character(1)` indicating how to collapse protein names for shared
  peptides. Default is `";"`.

- protColors:

  Either a `numeric(1)` or a named
  [`character()`](https://rdrr.io/r/base/character.html) of colour
  names. The numeric value indicates the protein colouring level to use.
  If 0 (default), all protein nodes are labelled in steelblue. For
  values \> 0, approximate string distances (see
  [`adist()`](https://rdrr.io/r/utils/adist.html)) between protein names
  are calculated and nodes of proteins that have names that differ will
  be coloured differently, with higher values leading to more colours.
  While no maximum to this value is defined in the code, it shouldn't be
  higher than the number of proteins. If a character is used, it should
  be a character of colour names named by protein identifiers. That
  vector should provide colours for at least all proteins in the
  adjacency matrix `m`, but more protein could be named. The latter is
  useful when generating a colour vector for all proteins in a dataset
  and use it for different adjacency matrix visualisations.

- pepColors:

  Either `NULL` (default) for no peptide colouring (white nodes) or a
  named [`character()`](https://rdrr.io/r/base/character.html) of colour
  names. It should be a character of colour names named by peptide
  identifiers. That vector should provide colours for at least all
  peptides in the adjacency matrix `m`, but more peptides could be
  named. The latter is useful when generating a colour vector for all
  peptides in a dataset and use it for different adjacency matrix
  visualisations.

- layout:

  A graph layout, as defined in the `ipgraph` package. Default is
  [`igraph::layout_as_bipartite()`](https://r.igraph.org/reference/layout_as_bipartite.html).

## Value

A peptide-by-protein sparce adjacency matrix (or class `dgCMatrix` as
defined in the `Matrix` package) or peptide/protein vector.

## Details

The `makeAdjacencyMatrix()` function creates a peptide-by-protein
adjacency matrix from a `character` or an instance of class
[`PSM()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md).

The character is formatted as
`x <- c("ProtA", "ProtB", "ProtA;ProtB", ...)`, as commonly encountered
in proteomics data spreadsheets. It defines that the first peptide is
mapped to protein "ProtA", the second one to protein "ProtB", the third
one to "ProtA" and "ProtB", and so on. The resulting matrix contains as
many rows as there are unique peptides and as many columns as there are
unique protein identifiers in `x`. The columns are named after the
protein identifiers and the peptide/protein vector names are used to
name the matrix rows (retaining only the unique names).

The `makePeptideProteinVector()` function does the opposite operation,
taking an adjacency matrix as input and returning a peptide/protein
vector. The matrix colnames are used to populate the vector and the
matrix rownames are used to name the vector elements.

Note that when creating an adjacency matrix from a PSM object, the
matrix is not necessarily binary, as multiple PSMs can match the same
peptide (sequence), such as for example precursors with different charge
states. A binary matrix can either be generated with the `binary`
argument (setting all non-0 values to 1) or by reducing the PSM object
accordingly (see example below).

It is also possible to generate adjacency matrices populated with
identification scores or probabilites by setting the "score" PSM
variable upon construction of the PSM object (see
[`PSM()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md)
for details). In case multiple PSMs occur, their respective scores get
summed.

The `plotAdjacencyMatrix()` function is useful to visualise small
adjacency matrices, such as those representing protein groups modelled
as connected components, as described and illustrated in
[`ConnectedComponents()`](https://rformassspectrometry.github.io/PSMatch/reference/ConnectedComponents.md).
The function generates a graph modelling the relation between proteins
(represented as squares) and peptides (represented as circes), which can
further be coloured (see the `protColors` and `pepColors` arguments).
The function invisibly returns the graph `igraph` object for additional
tuning and/or interactive visualisation using, for example
[`igraph::tkplot()`](https://r.igraph.org/reference/tkplot.html).

Such as illustrated in the examples below, each row/peptide is expected
to refer to protein groups or individual proteins (groups of size 1).
These have to be split accordingly.

## Author

Laurent Gatto

## Examples

``` r

## -----------------------
## From a character
## -----------------------

## Protein vector without names
prots <- c("ProtA", "ProtB", "ProtA;ProtB")
makeAdjacencyMatrix(prots)
#> 3 x 2 sparse Matrix of class "dgCMatrix"
#>   ProtA ProtB
#> 1     1     .
#> 2     .     1
#> 3     1     1

## Named protein vector
names(prots) <- c("pep1", "pep2", "pep3")
prots
#>          pep1          pep2          pep3 
#>       "ProtA"       "ProtB" "ProtA;ProtB" 
m <- makeAdjacencyMatrix(prots)
m
#> 3 x 2 sparse Matrix of class "dgCMatrix"
#>      ProtA ProtB
#> pep1     1     .
#> pep2     .     1
#> pep3     1     1

## Back to vector
vec <- makePeptideProteinVector(m)
vec
#>          pep1          pep2          pep3 
#>       "ProtA"       "ProtB" "ProtA;ProtB" 
identical(prots, vec)
#> [1] TRUE

## ----------------------------
## PSM object from a data.frame

## ----------------------------
## Case 1: Duplicate identifications

psmdf <- data.frame(psm = paste0("psm", 1:10),
                    peptide = paste0("pep", c(1, 1, 2, 2, 3, 4, 6, 7, 8, 8)),
                    protein = paste0("Prot", LETTERS[c(1, 1, 2, 2, 3, 4, 3, 5, 6, 6)]))
psmdf
#>      psm peptide protein
#> 1   psm1    pep1   ProtA
#> 2   psm2    pep1   ProtA
#> 3   psm3    pep2   ProtB
#> 4   psm4    pep2   ProtB
#> 5   psm5    pep3   ProtC
#> 6   psm6    pep4   ProtD
#> 7   psm7    pep6   ProtC
#> 8   psm8    pep7   ProtE
#> 9   psm9    pep8   ProtF
#> 10 psm10    pep8   ProtF
psm <- PSM(psmdf, peptide = "peptide", protein = "protein")
psm
#> PSM with 10 rows and 3 columns.
#> names(3): psm peptide protein
makeAdjacencyMatrix(psm)
#> 7 x 6 sparse Matrix of class "dgCMatrix"
#>      ProtA ProtB ProtC ProtD ProtE ProtF
#> pep1     2     .     .     .     .     .
#> pep2     .     2     .     .     .     .
#> pep3     .     .     1     .     .     .
#> pep4     .     .     .     1     .     .
#> pep6     .     .     1     .     .     .
#> pep7     .     .     .     .     1     .
#> pep8     .     .     .     .     .     2

## Reduce PSM object to peptides
rpsm <- reducePSMs(psm, k = psm$peptide)
rpsm
#> Reduced PSM with 7 rows and 3 columns.
#> names(3): psm peptide protein
makeAdjacencyMatrix(rpsm)
#> 7 x 6 sparse Matrix of class "dgCMatrix"
#>      ProtA ProtB ProtC ProtD ProtE ProtF
#> pep1     1     .     .     .     .     .
#> pep2     .     1     .     .     .     .
#> pep3     .     .     1     .     .     .
#> pep4     .     .     .     1     .     .
#> pep6     .     .     1     .     .     .
#> pep7     .     .     .     .     1     .
#> pep8     .     .     .     .     .     1

## Or set binary to TRUE
makeAdjacencyMatrix(psm, binary = TRUE)
#> 7 x 6 sparse Matrix of class "dgCMatrix"
#>      ProtA ProtB ProtC ProtD ProtE ProtF
#> pep1     1     .     .     .     .     .
#> pep2     .     1     .     .     .     .
#> pep3     .     .     1     .     .     .
#> pep4     .     .     .     1     .     .
#> pep6     .     .     1     .     .     .
#> pep7     .     .     .     .     1     .
#> pep8     .     .     .     .     .     1

## ----------------------------
## Case 2: Protein groups are separated by a semicolon
psmdf <- data.frame(psm = paste0("psm", 1:5),
                    peptide = paste0("pep", c(1, 2, 3, 4, 5)),
                    protein = c("ProtA", "ProtB;ProtD", "ProtA;ProtC",
                                "ProtC", "ProtA;ProtC;ProtD"))
psmdf
#>    psm peptide           protein
#> 1 psm1    pep1             ProtA
#> 2 psm2    pep2       ProtB;ProtD
#> 3 psm3    pep3       ProtA;ProtC
#> 4 psm4    pep4             ProtC
#> 5 psm5    pep5 ProtA;ProtC;ProtD
psm <- PSM(psmdf, peptide = "peptide", protein = "protein")
psm
#> PSM with 5 rows and 3 columns.
#> names(3): psm peptide protein
makeAdjacencyMatrix(psm, split = ";")
#> 5 x 4 sparse Matrix of class "dgCMatrix"
#>      ProtA ProtB ProtD ProtC
#> pep1     1     .     .     .
#> pep2     .     1     1     .
#> pep3     1     .     .     1
#> pep4     .     .     .     1
#> pep5     1     .     1     1

## ----------------------------
## PSM object from an mzid file
## ----------------------------

f <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache

psm <- PSM(f) |>
       filterPsmDecoy() |>
       filterPsmRank()
#> Removed 2896 decoy hits.
#> Removed 155 PSMs with rank > 1.
psm
#> PSM with 2751 rows and 35 columns.
#> names(35): sequence spectrumID ... subReplacementResidue subLocation
adj <- makeAdjacencyMatrix(psm)
dim(adj)
#> [1] 2357 1504
adj[1:10, 1:4]
#> 10 x 4 sparse Matrix of class "dgCMatrix"
#>                                   ECA2006 ECA1676 ECA3009 ECA1420
#> RQCRTDFLNYLR                            1       .       .       .
#> ESVALADQVTCVDWRNRKATKK                  .       1       .       .
#> QRMARTSDKQQSIRFLERLCGR                  .       .       2       .
#> DGGPAIYGHERVGRNAKKFKCLKFR               .       .       .       1
#> CIDRARHVEVQIFGDGKGRVVALGERDCSLQRR       .       .       .       .
#> VGRCRPIINYLASPGGTSER                    .       .       .       .
#> QRLDEHCVGVGQNALLLGR                     .       .       .       .
#> VDYQGKKVVIIGLGLTGLSCVDFFLARGVVPR        .       .       .       .
#> TLIIERCR                                .       .       .       .
#> QEVESCTTIEKIWVIDLKIK                    .       .       .       .

## Binary adjacency matrix
adj <- makeAdjacencyMatrix(psm, binary = TRUE)
adj[1:10, 1:4]
#> 10 x 4 sparse Matrix of class "dgCMatrix"
#>                                   ECA2006 ECA1676 ECA3009 ECA1420
#> RQCRTDFLNYLR                            1       .       .       .
#> ESVALADQVTCVDWRNRKATKK                  .       1       .       .
#> QRMARTSDKQQSIRFLERLCGR                  .       .       1       .
#> DGGPAIYGHERVGRNAKKFKCLKFR               .       .       .       1
#> CIDRARHVEVQIFGDGKGRVVALGERDCSLQRR       .       .       .       .
#> VGRCRPIINYLASPGGTSER                    .       .       .       .
#> QRLDEHCVGVGQNALLLGR                     .       .       .       .
#> VDYQGKKVVIIGLGLTGLSCVDFFLARGVVPR        .       .       .       .
#> TLIIERCR                                .       .       .       .
#> QEVESCTTIEKIWVIDLKIK                    .       .       .       .

## Peptides with rowSums > 1 match multiple proteins.
## Use filterPsmShared() to filter these out.
table(Matrix::rowSums(adj))
#> 
#>    1    2    3    4 
#> 2324   28    3    2 

## -----------------------------------------------
## Binary, non-binary and score adjacency matrices
## -----------------------------------------------

## -------------------------------------
## Case 1: no scores, 1 PSM per peptides
psmdf <- data.frame(spectrum = c("sp1", "sp2", "sp3", "sp4", "sp5",
                                 "sp6", "sp7", "sp8", "sp9", "sp10"),
                    sequence = c("NKAVRTYHEQ", "IYNHSQGFCA", "YHWRLPVSEF",
                                 "YEHNGFPLKD", "WAQFDVYNLS", "EDHINCTQWP",
                                 "WSMKVDYEQT", "GWTSKMRYPL", "PMAYIWEKLC",
                                 "HWAEYFNDVT"),
                    protein = c("ProtB", "ProtB", "ProtA", "ProtD", "ProtD",
                                "ProtG", "ProtF", "ProtE", "ProtC", "ProtF"),
                    decoy = rep(FALSE, 10),
                    rank = rep(1, 10),
                    score = c(0.082, 0.310, 0.133, 0.174, 0.944, 0.0261,
                              0.375, 0.741, 0.254, 0.058))
psmdf
#>    spectrum   sequence protein decoy rank  score
#> 1       sp1 NKAVRTYHEQ   ProtB FALSE    1 0.0820
#> 2       sp2 IYNHSQGFCA   ProtB FALSE    1 0.3100
#> 3       sp3 YHWRLPVSEF   ProtA FALSE    1 0.1330
#> 4       sp4 YEHNGFPLKD   ProtD FALSE    1 0.1740
#> 5       sp5 WAQFDVYNLS   ProtD FALSE    1 0.9440
#> 6       sp6 EDHINCTQWP   ProtG FALSE    1 0.0261
#> 7       sp7 WSMKVDYEQT   ProtF FALSE    1 0.3750
#> 8       sp8 GWTSKMRYPL   ProtE FALSE    1 0.7410
#> 9       sp9 PMAYIWEKLC   ProtC FALSE    1 0.2540
#> 10     sp10 HWAEYFNDVT   ProtF FALSE    1 0.0580

psm <- PSM(psmdf, spectrum = "spectrum", peptide = "sequence",
           protein = "protein", decoy = "decoy", rank = "rank")

## binary matrix
makeAdjacencyMatrix(psm)
#> 10 x 7 sparse Matrix of class "dgCMatrix"
#>            ProtB ProtA ProtD ProtG ProtF ProtE ProtC
#> NKAVRTYHEQ     1     .     .     .     .     .     .
#> IYNHSQGFCA     1     .     .     .     .     .     .
#> YHWRLPVSEF     .     1     .     .     .     .     .
#> YEHNGFPLKD     .     .     1     .     .     .     .
#> WAQFDVYNLS     .     .     1     .     .     .     .
#> EDHINCTQWP     .     .     .     1     .     .     .
#> WSMKVDYEQT     .     .     .     .     1     .     .
#> GWTSKMRYPL     .     .     .     .     .     1     .
#> PMAYIWEKLC     .     .     .     .     .     .     1
#> HWAEYFNDVT     .     .     .     .     1     .     .

## Case 2: sp1 and sp11 match the same peptide (NKAVRTYHEQ) as different PSMs
psmdf2 <- rbind(psmdf,
                data.frame(spectrum = "sp11",
                           sequence = psmdf$sequence[1],
                           protein = psmdf$protein[1],
                           decoy = FALSE,
                           rank = 1,
                           score = 0.011))
psmdf2
#>    spectrum   sequence protein decoy rank  score
#> 1       sp1 NKAVRTYHEQ   ProtB FALSE    1 0.0820
#> 2       sp2 IYNHSQGFCA   ProtB FALSE    1 0.3100
#> 3       sp3 YHWRLPVSEF   ProtA FALSE    1 0.1330
#> 4       sp4 YEHNGFPLKD   ProtD FALSE    1 0.1740
#> 5       sp5 WAQFDVYNLS   ProtD FALSE    1 0.9440
#> 6       sp6 EDHINCTQWP   ProtG FALSE    1 0.0261
#> 7       sp7 WSMKVDYEQT   ProtF FALSE    1 0.3750
#> 8       sp8 GWTSKMRYPL   ProtE FALSE    1 0.7410
#> 9       sp9 PMAYIWEKLC   ProtC FALSE    1 0.2540
#> 10     sp10 HWAEYFNDVT   ProtF FALSE    1 0.0580
#> 11     sp11 NKAVRTYHEQ   ProtB FALSE    1 0.0110
psm2 <- PSM(psmdf2, spectrum = "spectrum", peptide = "sequence",
            protein = "protein", decoy = "decoy", rank = "rank")

## Now NKAVRTYHEQ/ProtB counts 2 PSMs
makeAdjacencyMatrix(psm2)
#> 10 x 7 sparse Matrix of class "dgCMatrix"
#>            ProtB ProtA ProtD ProtG ProtF ProtE ProtC
#> NKAVRTYHEQ     2     .     .     .     .     .     .
#> IYNHSQGFCA     1     .     .     .     .     .     .
#> YHWRLPVSEF     .     1     .     .     .     .     .
#> YEHNGFPLKD     .     .     1     .     .     .     .
#> WAQFDVYNLS     .     .     1     .     .     .     .
#> EDHINCTQWP     .     .     .     1     .     .     .
#> WSMKVDYEQT     .     .     .     .     1     .     .
#> GWTSKMRYPL     .     .     .     .     .     1     .
#> PMAYIWEKLC     .     .     .     .     .     .     1
#> HWAEYFNDVT     .     .     .     .     1     .     .

## Force a binary matrix
makeAdjacencyMatrix(psm2, binary = TRUE)
#> 10 x 7 sparse Matrix of class "dgCMatrix"
#>            ProtB ProtA ProtD ProtG ProtF ProtE ProtC
#> NKAVRTYHEQ     1     .     .     .     .     .     .
#> IYNHSQGFCA     1     .     .     .     .     .     .
#> YHWRLPVSEF     .     1     .     .     .     .     .
#> YEHNGFPLKD     .     .     1     .     .     .     .
#> WAQFDVYNLS     .     .     1     .     .     .     .
#> EDHINCTQWP     .     .     .     1     .     .     .
#> WSMKVDYEQT     .     .     .     .     1     .     .
#> GWTSKMRYPL     .     .     .     .     .     1     .
#> PMAYIWEKLC     .     .     .     .     .     .     1
#> HWAEYFNDVT     .     .     .     .     1     .     .

## Case 3: Peptide (NKAVRTYHEQ) stems from multiple proteins (ProtB and
## ProtG). They are separated by a semicolon.
psmdf3 <- psmdf
psmdf3[psmdf3$sequence == "NKAVRTYHEQ","protein"] <- "ProtB;ProtG"
psmdf3
#>    spectrum   sequence     protein decoy rank  score
#> 1       sp1 NKAVRTYHEQ ProtB;ProtG FALSE    1 0.0820
#> 2       sp2 IYNHSQGFCA       ProtB FALSE    1 0.3100
#> 3       sp3 YHWRLPVSEF       ProtA FALSE    1 0.1330
#> 4       sp4 YEHNGFPLKD       ProtD FALSE    1 0.1740
#> 5       sp5 WAQFDVYNLS       ProtD FALSE    1 0.9440
#> 6       sp6 EDHINCTQWP       ProtG FALSE    1 0.0261
#> 7       sp7 WSMKVDYEQT       ProtF FALSE    1 0.3750
#> 8       sp8 GWTSKMRYPL       ProtE FALSE    1 0.7410
#> 9       sp9 PMAYIWEKLC       ProtC FALSE    1 0.2540
#> 10     sp10 HWAEYFNDVT       ProtF FALSE    1 0.0580
psm3 <- PSM(psmdf3, spectrum = "spectrum", peptide = "sequence",
            protein = "protein", decoy = "decoy", rank = "rank")

## Now ProtB & ProtG count 2 PSMs each: NKAVRTYHEQ and IYNHSQGFCA &
## EDHINCTQWP respectively
makeAdjacencyMatrix(psm3, split = ";")
#> 10 x 7 sparse Matrix of class "dgCMatrix"
#>            ProtB ProtG ProtA ProtD ProtF ProtE ProtC
#> NKAVRTYHEQ     1     1     .     .     .     .     .
#> IYNHSQGFCA     1     .     .     .     .     .     .
#> YHWRLPVSEF     .     .     1     .     .     .     .
#> YEHNGFPLKD     .     .     .     1     .     .     .
#> WAQFDVYNLS     .     .     .     1     .     .     .
#> EDHINCTQWP     .     1     .     .     .     .     .
#> WSMKVDYEQT     .     .     .     .     1     .     .
#> GWTSKMRYPL     .     .     .     .     .     1     .
#> PMAYIWEKLC     .     .     .     .     .     .     1
#> HWAEYFNDVT     .     .     .     .     1     .     .

## --------------------------------
## Case 4: set the score PSM values
psmVariables(psm) ## no score defined
#>   spectrum    peptide    protein      decoy       rank      score        fdr 
#> "spectrum" "sequence"  "protein"    "decoy"     "rank"         NA         NA 
psm4 <- PSM(psm, spectrum = "spectrum", peptide = "sequence",
            protein = "protein", decoy = "decoy", rank = "rank",
            score = "score")
psmVariables(psm4) ## score defined
#>   spectrum    peptide    protein      decoy       rank      score        fdr 
#> "spectrum" "sequence"  "protein"    "decoy"     "rank"    "score"         NA 

## adjacency matrix with scores
makeAdjacencyMatrix(psm4)
#> 10 x 7 sparse Matrix of class "dgCMatrix"
#>            ProtB ProtA ProtD  ProtG ProtF ProtE ProtC
#> NKAVRTYHEQ 0.082 .     .     .      .     .     .    
#> IYNHSQGFCA 0.310 .     .     .      .     .     .    
#> YHWRLPVSEF .     0.133 .     .      .     .     .    
#> YEHNGFPLKD .     .     0.174 .      .     .     .    
#> WAQFDVYNLS .     .     0.944 .      .     .     .    
#> EDHINCTQWP .     .     .     0.0261 .     .     .    
#> WSMKVDYEQT .     .     .     .      0.375 .     .    
#> GWTSKMRYPL .     .     .     .      .     0.741 .    
#> PMAYIWEKLC .     .     .     .      .     .     0.254
#> HWAEYFNDVT .     .     .     .      0.058 .     .    

## Force a binary matrix
makeAdjacencyMatrix(psm4, binary = TRUE)
#> 10 x 7 sparse Matrix of class "dgCMatrix"
#>            ProtB ProtA ProtD ProtG ProtF ProtE ProtC
#> NKAVRTYHEQ     1     .     .     .     .     .     .
#> IYNHSQGFCA     1     .     .     .     .     .     .
#> YHWRLPVSEF     .     1     .     .     .     .     .
#> YEHNGFPLKD     .     .     1     .     .     .     .
#> WAQFDVYNLS     .     .     1     .     .     .     .
#> EDHINCTQWP     .     .     .     1     .     .     .
#> WSMKVDYEQT     .     .     .     .     1     .     .
#> GWTSKMRYPL     .     .     .     .     .     1     .
#> PMAYIWEKLC     .     .     .     .     .     .     1
#> HWAEYFNDVT     .     .     .     .     1     .     .

## ---------------------------------
## Case 5: scores with multiple PSMs

psm5 <- PSM(psm2, spectrum = "spectrum", peptide = "sequence",
            protein = "protein", decoy = "decoy", rank = "rank",
            score = "score")

## Now NKAVRTYHEQ/ProtB has a summed score of 0.093 computed as
## 0.082 (from sp1) + 0.011 (from sp11)
makeAdjacencyMatrix(psm5)
#> 10 x 7 sparse Matrix of class "dgCMatrix"
#>            ProtB ProtA ProtD  ProtG ProtF ProtE ProtC
#> NKAVRTYHEQ 0.093 .     .     .      .     .     .    
#> IYNHSQGFCA 0.310 .     .     .      .     .     .    
#> YHWRLPVSEF .     0.133 .     .      .     .     .    
#> YEHNGFPLKD .     .     0.174 .      .     .     .    
#> WAQFDVYNLS .     .     0.944 .      .     .     .    
#> EDHINCTQWP .     .     .     0.0261 .     .     .    
#> WSMKVDYEQT .     .     .     .      0.375 .     .    
#> GWTSKMRYPL .     .     .     .      .     0.741 .    
#> PMAYIWEKLC .     .     .     .      .     .     0.254
#> HWAEYFNDVT .     .     .     .      0.058 .     .    
```

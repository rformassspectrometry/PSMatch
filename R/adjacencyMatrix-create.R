##' @title Convert to/from an adjacency matrix.
##'
##' @description
##'
##' There are two ways that peptide/protein matches are commonly
##' stored: either as a vector or an adjacency matrix. The functions
##' described below convert between these two format.
##'
##' @details
##'
##' The [makeAdjacencyMatrix()] function creates a peptide-by-protein adjacency
##' matrix from a `character` or an instance of class [PSM()].
##'
##' The character is formatted as `x <- c("ProtA", "ProtB", "ProtA;ProtB",
##' ...)`, as commonly encoutered in proteomics data spreadsheets. It defines
##' that the first peptide is mapped to protein "ProtA", the second one to
##' protein "ProtB", the third one to "ProtA" and "ProtB", and so on. The
##' resulting matrix contains `length(x)` rows and as many columns as there are
##' unique protein idenifiers in `x`. The columns are named after the protein
##' identifiers and the peptide/protein vector names are used to name to matrix
##' rows (even if these aren't unique).
##'
##' The [makePeptideProteinVector()] function does the opposite operation,
##' taking an adjacency matrix as input and retruning a peptide/protein
##' vector. The matrix colnames are used to populate the vector and the matrix
##' rownames are used to name the vector elements.
##'
##' Note that when creating an adjacency matrix from a PSM object, the matrix is
##' not necessarily binary, as multiple PSMs can match the same peptide
##' (sequence), such as for example precursors with different charge states. A
##' binary matrix can either be generated with the `binary` argument (setting
##' all non-0 values to 1) or by reducing the PSM object accordingly (see
##' example below).
##'
##' It is also possible to generate adjacency matrices populated with
##' identification scores or probabilites by setting the "score" PSM variable
##' upon construction of the PSM object (see [PSM()] for details). In case
##' multiple PSMs occur, their respective scores get summed.
##'
##' The `plotAdjacencyMatrix()` function is useful to visualise small adjacency
##' matrices, such as those representing protein groups modelled as connected
##' components, as described and illustrated in [ConnectedComponents()]. The
##' function generates a graph modelling the relation between proteins
##' (represented as squares) and peptides (represented as circes), which can
##' further be coloured (see the `protColors` and `pepColors` arguments). The
##' function invisibly returns the graph `igraph` object for additional tuning
##' and/or interactive visualisation using, for example [igraph::tkplot()].
##'
##' There exists some important differences in the creation of an adjacency
##' matrix from a PSM object or a vector, other than the input variable itself:
##'
##' - In a `PSM` object, each row (PSM) refers to an *individual* proteins;
##'   rows/PSMs never refer to a protein group. There is thus no need for a
##'   `split` argument, which is used exclusively when creating a matrix from a
##'   character.
##'
##' - Conversely, when using protein vectors, such as those illustrated in the
##'   example below or retrieved from tabular quantitative proteomics data, each
##'   row/peptide is expected to refer to protein groups or individual proteins
##'   (groups of size 1). These have to be split accordingly.
##'
##' @param x Either an instance of class `PSM` or a `character`. See
##'     example below for details.
##'
##' @param split `character(1)` defining how to split the string of
##'     protein identifiers (using [strsplit()]). Default is ";". If
##'     `NULL`, splitting is ignored.
##'
##' @param peptide `character(1)` indicating the name of the variable
##'     that defines peptides in the `PSM` object. Default is the
##'     `peptide` PSM variable as defined in [psmVariables()].
##'
##' @param protein `character(1)` indicating the name of the variable
##'     that defines proteins in the `PSM` object. Default is the
##'     `peptide` PSM variable as defined in [psmVariables()].
##'
##' @param score `character(1)` indicating the name of the variable
##'     that defines PSM scores in the `PSM` object. Default is the
##'     `score` PSM variable as defined in [psmVariables()]. Ignored
##'     when `NA` (which is the default value unless set by the user
##'     when constructing the `PSM` object).
##'
##' @param binary `logical(1)` indicates if the adjacency matrix
##'     should be strictly binary. In such a case, PSMs matching the
##'     same peptide but from different precursors (for example charge
##'     2 and 3) or carrying different PTMs, are counted only
##'     once. Default if `FALSE`. This also overrides any `score` that
##'     would be set.
##'
##' @return A peptide-by-protein sparce adjacency matrix (or class
##'     `dgCMatrix` as defined in the `Matrix` package) or
##'     peptide/protein vector.
##'
##' @author Laurent Gatto
##'
##' @name adjacencyMatrix
##'
##' @export
##'
##' @importFrom Matrix Matrix sparseMatrix
##'
##' @examples
##'
##' ## -----------------------
##' ## From a character
##' ## -----------------------
##'
##' ## Protein vector without names
##' prots <- c("ProtA", "ProtB", "ProtA;ProtB")
##' makeAdjacencyMatrix(prots)
##'
##' ## Named protein vector
##' names(prots) <- c("pep1", "pep2", "pep3")
##' prots
##' m <- makeAdjacencyMatrix(prots)
##' m
##'
##' ## Back to vector
##' vec <- makePeptideProteinVector(m)
##' vec
##' identical(prots, vec)
##'
##' ## ----------------------------
##' ## PSM object from a data.frame
##' ## ----------------------------
##'
##' psmdf <- data.frame(psm = paste0("psm", 1:10),
##'                     peptide = paste0("pep", c(1, 1, 2, 2, 3, 4, 6, 7, 8, 8)),
##'                     protein = paste0("Prot", LETTERS[c(1, 1, 2, 2, 3, 4, 3, 5, 6, 6)]))
##' psmdf
##' psm <- PSM(psmdf, peptide = "peptide", protein = "protein")
##' psm
##' makeAdjacencyMatrix(psm)
##'
##' ## Reduce PSM object to peptides
##' rpsm <- reducePSMs(psm, k = psm$peptide)
##' rpsm
##' makeAdjacencyMatrix(rpsm)
##'
##' ## Or set binary to TRUE
##' makeAdjacencyMatrix(psm, binary = TRUE)
##'
##' ## ----------------------------
##' ## PSM object from an mzid file
##' ## ----------------------------
##'
##' f <- msdata::ident(full.names = TRUE, pattern = "TMT")
##' psm <- PSM(f) |>
##'        filterPsmDecoy() |>
##'        filterPsmRank()
##' psm
##' adj <- makeAdjacencyMatrix(psm)
##' dim(adj)
##' adj[1:10, 1:4]
##'
##' ## Binary adjacency matrix
##' adj <- makeAdjacencyMatrix(psm, binary = TRUE)
##' adj[1:10, 1:4]
##'
##' ## Peptides with rowSums > 1 match multiple proteins.
##' ## Use filterPsmShared() to filter these out.
##' table(Matrix::rowSums(adj))
##'
##' ## -----------------------------------------------
##' ## Binary, non-binary and score adjacency matrices
##' ## -----------------------------------------------
##'
##' ## -------------------------------------
##' ## Case 1: no scores, 1 PSM per peptides
##' psmdf <- data.frame(spectrum = c("sp1", "sp2", "sp3", "sp4", "sp5",
##'                                  "sp6", "sp7", "sp8", "sp9", "sp10"),
##'                     sequence = c("NKAVRTYHEQ", "IYNHSQGFCA", "YHWRLPVSEF",
##'                                  "YEHNGFPLKD", "WAQFDVYNLS", "EDHINCTQWP",
##'                                  "WSMKVDYEQT", "GWTSKMRYPL", "PMAYIWEKLC",
##'                                  "HWAEYFNDVT"),
##'                     protein = c("ProtB", "ProtB", "ProtA", "ProtD", "ProtD",
##'                                 "ProtG", "ProtF", "ProtE", "ProtC", "ProtF"),
##'                     decoy = rep(FALSE, 10),
##'                     rank = rep(1, 10),
##'                     score = c(0.082, 0.310, 0.133, 0.174, 0.944, 0.0261,
##'                               0.375, 0.741, 0.254, 0.058))
##' psmdf
##'
##' psm <- PSM(psmdf, spectrum = "spectrum", peptide = "sequence",
##'            protein = "protein", decoy = "decoy", rank = "rank")
##'
##' ## binary matrix
##' makeAdjacencyMatrix(psm)
##'
##  ## --------------------------------------------------------
##' ## Case 2: sp1 and sp11 match the same peptide (NKAVRTYHEQ)
##' psmdf2 <- rbind(psmdf,
##'                 data.frame(spectrum = "sp11",
##'                            sequence = psmdf$sequence[1],
##'                            protein = psmdf$protein[1],
##'                            decoy = FALSE,
##'                            rank = 1,
##'                            score = 0.011))
##' psmdf2
##' psm2 <- PSM(psmdf2, spectrum = "spectrum", peptide = "sequence",
##'             protein = "protein", decoy = "decoy", rank = "rank")
##'
##' ## Now NKAVRTYHEQ/ProtB counts 2 PSMs
##' makeAdjacencyMatrix(psm2)
##'
##' ## Force a binary matrix
##' makeAdjacencyMatrix(psm2, binary = TRUE)
##'
##' ## --------------------------------
##' ## Case 3: set the score PSM values
##' psmVariables(psm) ## no score defined
##' psm3 <- PSM(psm, spectrum = "spectrum", peptide = "sequence",
##'             protein = "protein", decoy = "decoy", rank = "rank",
##'             score = "score")
##' psmVariables(psm3) ## score defined
##'
##' ## adjacency matrix with scores
##' makeAdjacencyMatrix(psm3)
##'
##' ## Force a binary matrix
##' makeAdjacencyMatrix(psm3, binary = TRUE)
##'
##' ## ---------------------------------
##' ## Case 4: scores with multiple PSMs
##'
##' psm4 <- PSM(psm2, spectrum = "spectrum", peptide = "sequence",
##'             protein = "protein", decoy = "decoy", rank = "rank",
##'             score = "score")
##'
##' ## Now NKAVRTYHEQ/ProtB has a summed score of 0.093 computed as
##' ## 0.082 (from sp1) + 0.011 (from sp11)
##' makeAdjacencyMatrix(psm4)
makeAdjacencyMatrix <- function(x, split = ";",
                                peptide = psmVariables(x)["peptide"],
                                protein = psmVariables(x)["protein"],
                                score = psmVariables(x)["score"],
                                binary = FALSE) {
    if (inherits(x, "PSM")) {
        adj <- .makeSparseAdjacencyMatrixFromPSM(x, peptide, protein, score, split)
    } else if (is.character(x)) {
        adj <- .makeSparseAdjacencyMatrixFromChar(x, split)
    } else stop("'x' must be a character or a PSM object.")
    if (binary)
        adj[abs(adj) > 0] <- 1
    return(adj)
}

.makeSparseAdjacencyMatrixFromChar <- function(x, split = ";") {
    if (is.null(split)) col_list <- x
    else col_list <- strsplit(x, split)
    if (is.null(names(col_list))) {
            col_list_names <- seq_along(col_list)
    } else col_list_names <- names(col_list)
    row_names <- unique(col_list_names)
    col_names <- unique(unlist(col_list))
    i <- rep(col_list_names, lengths(col_list))
    i <- match(i, row_names)
    j <- unname(unlist(lapply(col_list, match, col_names)))
    sparseMatrix(i, j, x = 1, dimnames = list(row_names, col_names))
}

.makeSparseAdjacencyMatrixFromPSM <- function(x, peptide, protein, score) {
    if (is.na(peptide) | is.na(protein))
        stop("Please define the 'protein' and 'peptide' PSM variables.")
    if (!protein %in% names(x) | !peptide %in% names(x))
        stop("PSM variables 'protein' and 'peptide' must be defined.")
    row_names <- unique(x[[peptide]])
    col_names <- unique(x[[protein]])
    i <- match(x[[peptide]], row_names)
    j <- match(x[[protein]], col_names)
    adj_values <- 1
     if (!is.na(score))
        adj_values <- x[[score]]
    sparseMatrix(i, j, x = adj_values,
                 dimnames = list(row_names, col_names))
}


##' @param m A peptide-by-protein adjacency matrix.
##'
##' @param collapse `character(1)` indicating how to collapse protein
##'     names for shared peptides. Default is `";"`.
##'
##' @name adjacencyMatrix
##'
##' @export
makePeptideProteinVector <- function(m, collapse = ";") {
    stopifnot(is.matrix(m) | inherits(m, "Matrix"))
    vec <- rep(NA_character_, nrow(m))
    for (i in seq_len(nrow(m)))
        vec[i] <- paste(names(which(m[i, ] != 0)), collapse = collapse)
    names(vec) <- rownames(m)
    vec
}
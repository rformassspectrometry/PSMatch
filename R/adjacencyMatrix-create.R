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
##' The [makeAdjacencyMatrix()] function creates a peptide-by-protein
##' adjacency matrix from a `character` or an instance of class
##' [PSM()].
##'
##' The character is formatted as `x <- c("ProtA", "ProtB",
##' "ProtA;ProtB", ...)`, as commonly encoutered in proteomics data
##' spreadsheets. It defines that the first peptide is mapped to
##' protein "ProtA", the second one to protein "ProtB", the third one
##' to "ProtA" and "ProtB", and so on. The resulting matrix contain
##' `length(x)` rows an as many columns as there are unique protein
##' idenifiers in `x`. The columns are named after the protein
##' idenifiers and the peptide/protein vector namesa are used to name
##' to matrix rows (even if these aren't unique).
##'
##' The [makePeptideProteinVector()] function does the opposite
##' operation, taking an adjacency matrix as input and retruning a
##' peptide/protein vector. The matrix colnames are used to populate
##' the vector and the matrix rownames are used to name the vector
##' elements.
##'
##' Note that when creating an adjacency matrix from PSMs, the matrix
##' is not necessarily binary, as multiple PSMs can match the same
##' peptide (sequence), such as for example precursors with different
##' charge states. A binary matrix can either be generated with the
##' `binary` argument (setting all non-0 values to 1) or by reducing
##' the PSM object accordingly (see example below).
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
##'     `peptide` PSM variable as defined in [psmVariables()].Default
##'     is `DatanbaseAccess`.
##'
##' @param binary `logical(1)` indicates if the adjacency matrix
##'     should be strictly binary. In such case, PSMs matching the
##'     same peptide but from different precursors (for example charge
##'     2 and 3) or carrying different PTMs, are counted only
##'     once. Default if `FALSE`.
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
##' ## ---------------------------
##' ## PSM object from a dat.frame
##' ## ---------------------------
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
##' table(rowSums(adj))
makeAdjacencyMatrix <- function(x, split = ";",
                                peptide = psmVariables(x)["peptide"],
                                protein = psmVariables(x)["protein"],
                                binary = FALSE) {
    if (inherits(x, "PSM")) {
        adj <- .makeSparseAdjacencyMatrixFromPSM(x, peptide, protein)
    } else if (is.character(x)) {
        adj <- .makeSparseAdjacencyMatrixFromChar(x, split)
    } else stop("'x' must be a character or a PSM object.")
    if (binary)
        adj[adj > 1] <- 1
    return(adj)
}

.makeSparseAdjacencyMatrixFromChar <- function(x, split = ";") {
    if (is.null(split)) col_list <- x
    else col_list <- strsplit(x, split)
    if (is.null(names(col_list))) {
        row_names <- seq_along(col_list)
    } else row_names <- names(col_list)
    col_names <- unique(unlist(col_list))
    i <- rep(seq_along(col_list), lengths(col_list))
    j <- unname(unlist(sapply(col_list, match, col_names)))
    sparseMatrix(i, j, x = 1, dimnames = list(row_names, col_names))
}

.makeSparseAdjacencyMatrixFromPSM <- function(x, peptide, protein) {
    stopifnot(peptide %in% names(x))
    stopifnot(protein %in% names(x))
    row_names <- unique(x[[peptide]])
    col_names <- unique(x[[protein]])
    i <- match(x[[peptide]], row_names)
    j <- match(x[[protein]], col_names)
    sparseMatrix(i, j, x = 1, dimnames = list(row_names, col_names))
}


##' @param m An adjacency matrix.
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

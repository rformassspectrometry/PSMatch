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
##' @return A peptide-by-protein adjacency matrix or peptide/protein
##'     vector.
##'
##' @author Laurent Gatto
##'
##' @name adjacencyMatrix
##'
##' @export
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
##' ## -----------------------
##' ## From a PSM object
##' ## -----------------------
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
    if (inherits(x, "PSM"))
        return(.makeAdjacencyMatrixFromPSM(x, peptide, protein, binary))
    if (is.character(x))
        return(.makeAdjacencyMatrixFromChar(x, split))
    stop("'x' must be a character or a PSM object.")
}

.makeAdjacencyMatrixFromChar <- function(x, split = ";") {
    n <- length(x)
    if (is.null(split)) {
        col_list <- x
        m <- length(cnames <- unique(x))
    } else {
        col_list <- strsplit(x, split)
        m <- length(cnames <- unique(unlist(col_list)))
    }
    adj <- matrix(0, nrow = n, ncol = m,
                  dimnames = list(names(x), cnames))
    for (i in seq_along(col_list)) {
        adj[i, col_list[[i]]] <- 1
    }
    adj
}

.makeAdjacencyMatrixFromPSM <- function(x, peptide, protein, binary) {
    n <- length(nx <- unique(x[[peptide]]))
    m <- length(mx <- unique(x[[protein]]))
    adj <- matrix(0, nrow = n, ncol = m,
                  dimnames = list(nx, mx))
    for (k in x[[peptide]]) {
        i <- which(x[[peptide]] %in% k)
        adj[k, x[[protein]][i]] <- adj[k, x[[protein]][i]] + 1
    }
    if (binary)
        adj[adj > 1] <- 1
    adj
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
    stopifnot(is.matrix(m))
    vec <- rep(NA_character_, nrow(m))
    for (i in seq_len(nrow(m)))
        vec[i] <- paste(names(which(m[i, ] != 0)), collapse = collapse)
    names(vec) <- rownames(m)
    vec
}

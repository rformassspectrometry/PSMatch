##' @title Create an adjacency matrix.
##'
##' @description
##'
##' This function created a peptide-by-protein adjacency matrix from a
##' `character` or an instance of class `PSM`.
##'
##' The character is formatted as `x <- c("ProtA", "ProtB",
##' "ProtA;ProtB", ...)`, as commonly encoutered in proteomics data
##' spreadsheets. It defines that the first peptide is mapped to
##' protein "ProtA", the second one to protein "ProtB", the third one
##' to "ProtA" and "ProtB", and so on. The resulting matrix contain
##' `length(x)` rows an as many columns as there are unique protein
##' idenifiers in `x`. The column are always named after the protein
##' idenifiers. If the protein identifier vector is named and the
##' names are unique, these are then used to name to matrix rows.
##'
##' @param x Either an instance of class `PSM` or a `character`. See
##'     example below for details.
##'
##' @param split A `character(1)` that defines how to split the string
##'     of protein identifiers.
##'
##' @return A peptide-by-protein adjacency `matrix`.
##'
##' @author Laurent Gatto
##'
##' @export
##'
##' @examples
##' ## Protein vector without names
##' prots <- c("ProtA", "ProtB", "ProtA;ProtB")
##' makeAdjacencyMatrix(prots)
##'
##' ## Named protein vector
##' names(prots) <- c("pep1", "pep2", "pep3")
##' makeAdjacencyMatrix(prots)
makeAdjacencyMatrix <- function(x, split = ";",
                                peptides = "sequence",
                                proteins = "DatabaseAccess",
                                sparse = FALSE) {
    if (inherits(x, "PSM"))
        return(.makeAdjacencyMatrixFromPSM(x, peptides, proteins, sparse))
    if (is.character(x) & !sparse)
        return(.makeAdjacencyMatrixFromChar(x, split))
    if (is.character(x) & sparse)
        return(.makeSparseAdjacencyMatrixFromChar(x, split))
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
    adj <- matrix(0, nrow = n, ncol = m)
    colnames(adj) <- cnames
    if (!is.null(names(x)) & !anyDuplicated(names(x))) {
        rownames(adj) <- names(x)
    }
    for (i in seq_along(col_list)) {
        adj[i, col_list[[i]]] <- 1
    }
    adj
}

.makeAdjacencyMatrixFromPSM <- function(x, peptides, proteins, sparse) {
    vec <- x[[proteins]]
    names(vec) <- x[[peptides]]
    ## make names unique?
    makeAdjacencyMatrix(vec, split = NULL, sparse)
}

.makeSparseAdjacencyMatrixFromChar <- function(x, split) {
    stop("Not yet implemented")
}

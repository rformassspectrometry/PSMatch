##' @title Create an adjacency matrix.
##'
##' @description
##'
##' This function created a peptide-by-protein adjacency matrix from a
##' `character` or an instance of class `PSM`.
##'
##' @details
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
##' The adjacency matrix produced from a `PSM` object doesn't
##' represent a peptide-by-protein adjacency matrix, given that the
##' observations of a `PSM` object represent peptide-spectrum matches,
##' rather than peptides. It is possible to set the `unique` parameter
##' to return unique PSM occurences. Note however that this implicitly
##' assumes that the same peptides where matched to the same
##' protein(s), which is not explicitly verified.
##'
##' @param x Either an instance of class `PSM` or a `character`. See
##'     example below for details.
##'
##' @param split `character(1)` defining how to split the string of
##'     protein identifiers (using [strsplit()]). Default is ";". If
##'     `NULL`, splitting is ignored.
##'
##' @param peptides `character(1)` indicating the name of the variable
##'     that defines peptides in the `PSM` object. Default is
##'     `sequence`.
##'
##' @param proteins `character(1)` indicating the name of the variable
##'     that defines proteins in the `PSM` object. Default is
##'     `DatanbaseAccess`.
##'
##' @param unique `logical(1)` defining if all peptides should should
##'     be considered in the contingency matrix, or should duplicates
##'     be ingored (and only first occurences taken into account). The
##'     default is `FALSE` and all peptides names are made unique by
##'     appending sequence numbers to duplicates (see
##'     [make.unique()]).
##'
##' @return A peptide-by-protein adjacency matrix.
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
##'
##' ## From a PSM object
##' f <- msdata::ident(full.names = TRUE, pattern = "TMT")
##' psm <- filterPSMs(PSM(f))
##' psm
##' adj <- makeAdjacencyMatrix(psm)
##' dim(adj)
##' adj[1:10, 1:4]
##'
##' ## Drop duplicated peptides
##' adj <- makeAdjacencyMatrix(psm, unique = TRUE)
##' dim(adj)
##' adj[1:10, 1:4]
makeAdjacencyMatrix <- function(x, split = ";",
                                peptides = "sequence",
                                proteins = "DatabaseAccess",
                                unique = FALSE) {
    if (inherits(x, "PSM"))
        return(.makeAdjacencyMatrixFromPSM(x, peptides, proteins, unique))
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

.makeAdjacencyMatrixFromPSM <- function(x,
                                        peptides, proteins,
                                        unique) {
    vec <- x[[proteins]]
    names(vec) <- make.unique(x[[peptides]])
    adj <- makeAdjacencyMatrix(vec, split = NULL)
    if (unique)
        adj <- adj[!duplicated(x[[peptides]]), ]
    adj
}

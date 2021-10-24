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
##' library(magrittr)
##' psm <- PSM(f) %>%
##'        filterPsmDecoy() %>%
##'        filterPsmRank()
##' psm
##' adj <- makeAdjacencyMatrix(psm)
##' dim(adj)
##' adj[1:10, 1:4]
##' ## Peptides with rowSums > 1 match multiple proteins.
##' ## Use filterPsmNonProteotypic() to filter these out.
##' table(rowSums(adj))
makeAdjacencyMatrix <- function(x, split = ";",
                                peptide = psmVariables(x)["peptide"],
                                protein = psmVariables(x)["protein"]) {
    if (inherits(x, "PSM"))
        return(.makeAdjacencyMatrixFromPSM(x, peptide, protein))
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

.makeAdjacencyMatrixFromPSM <- function(x, peptide, protein) {
    n <- length(nx <- unique(x[[peptide]]))
    m <- length(mx <- unique(x[[protein]]))
    adj <- matrix(0, nrow = n, ncol = m,
                  dimnames = list(rownames = nx,
                                  colnames = mx))
    for (k in nx) {
        i <- which(x[[peptide]] %in% k)
        adj[k, x[[protein]][i]] <- 1
    }
    adj
}

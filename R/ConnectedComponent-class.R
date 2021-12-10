##' @title Connected components
##'
##' @aliases ConnectedComponents ConnectedComponents-class
##'
##' @name ConnectedComponents
##'
##' @description
##'
##' Connected components are a useful representation when exploring
##' identification data. They represent the relation between proteins
##' (the connected components) and how they form groups of proteins
##' defined by shared peptides.
##'
##' Connected components are stored as `Connectedcomponents` that can
##' be generated using the `ConnectedComponents()` function.
##'
##' @examples
##'
##' ## --------------------------------
##' ## From an adjacency matrix
##' ## --------------------------------
##' library(Matrix)
##' adj <- sparseMatrix(i = c(1, 2, 3, 3), j = c(1, 2, 3, 4), x = 1,
##'        dimnames = list(paste0("Pep", 1:3),
##'                        paste0("Prot", 1:4)))
##' adj
##' cc <- ConnectedComponents(adj)
##' cc
##'
##' ## --------------------------------
##' ## From an PSM object
##' ## --------------------------------
##' f <- msdata::ident(full.names = TRUE, pattern = "TMT")
##' f
##'
##' psm <- PSM(f) |>
##'        filterPsmDecoy() |>
##'        filterPsmRank()
##'
##' cc <- ConnectedComponents(psm)
##' cc
NULL

setClass("ConnectedComponents",
         slots = c(adjMatrix = "Matrix",
                   ccMatrix = "Matrix",
                   ccList = "List",
                   ccPeptides = "List"))

##' @importFrom methods new
##'
##' @importFrom Matrix t tcrossprod
##'
##' @rdname ConnectedComponents
##'
##' @export
##'
##' @param object An adjacency matrix class `Matrix` or an instance of
##'     class `PSM`.
ConnectedComponents <- function(object) {
    if (is(object, "PSM")) {
        adj <- makeAdjacencyMatrix(object)
    } else if (is(object, "Matrix")) {
        adj <- object
    } else stop("'object' must be of class 'PSM' or 'Matrix.")
    getCCpeptides <- function(cc, adj) {
        res <- adj[, cc, drop = FALSE]
        res <- res[rowSums(res) > 0, , drop = FALSE]
        res
    }
    cc <- Matrix::tcrossprod(t(adj))
    n <- ncol(cc)
    cc_pep <- cc_list <- vector("list", length = n)
    for (i in seq_len(n-1)) {
        j <- i:n
        k <- which(cc[i, j] != 0)
        cc_list[[i]] <- colnames(adj)[j][k]
        cc_pep[[i]] <- rownames(getCCpeptides(cc_list[[i]], adj))
    }
    sel <- lengths(cc_list) > 1
    new("ConnectedComponents",
        adjMatrix = adj,
        ccMatrix = cc,
        ccList = List(cc_list[sel]),
        ccPeptides = List(cc_pep[sel]))
}


setMethod("show", "ConnectedComponents",
          function(object) {
              cat(sprintf("An instance of class %s", class(object)), "\n")
              cat(" Number of componenents ", nrow(object@ccMatrix), "\n")
              cat(" CC with size > 1:\n")
              tab <- table(lengths(object@ccList))
              msg <- strwrap(paste(paste0(tab, "(", names(tab), ")"),
                                   collapse = " "))
              message(paste(" ", msg, collapse = "\n"))
          })

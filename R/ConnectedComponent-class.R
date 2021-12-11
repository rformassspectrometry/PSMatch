##' @title Connected components
##'
##' @aliases ConnectedComponents ConnectedComponents-class
##'
##' @name ConnectedComponents
##'
##' @aliases ConnectedComponents-class ConnectedComponents length,ConnectedComponents lengths,ConnectedComponents adjacencyMatrix,ConnectedComponents ccMatrix ccList ccPeptides show,ConnectedComponents connectedComponents
##'
##' @description
##'
##' Connected components are a useful representation when exploring
##' identification data. They represent the relation between proteins
##' (the connected components) and how they form groups of proteins as
##' defined by shared peptides.
##'
##' Connected components are stored as `ConnectedComponents` objects
##' that can be generated using the `ConnectedComponents()`
##' function.
##'
##' @slot adjMatrix The sparse adjacency matrix (class `Matrix`) of
##'     dimension *p* peptides by *m* proteins that was used to
##'     generate the object.
##'
##' @slot ccMatrix The sparse connected components matrix (class
##'     `Matrix`) of dimension *m* by *m* proteins.
##'
##' @slot ccList A `List` containing the protein names composing the
##'     respective connected composing.
##'
##' @slot ccPeptides A `List` containing the peptide names found in
##'     the proteins composing the respective connected composing.
##'
##' @section Creating and manipulating objects:
##'
##' - Instances of the class are created with the
##'   `ConnectedComponent()` constructor from a [PSM()] object or
##'   directly from a sparse adjacency matrix of class `Matrix`. Note
##'   that if using the latter, the rows and columns must be named.
##'
##' - The sparse peptide-by-protein adjacency matrix is stored in the
##'   `ConnectedComponent` instance and can be accessed with the
##'   `adjacencyMatrix()` function.
##'
##' - The protein-by-protein connected components sparse matrix of
##'   object `x` can be accessed with the `ccMatrix(x)` function.
##'
##' - The `ccList(x)` function returns a `List` with the proteins
##'   connected components of object `x`.
##'
##' - The `ccPeptides(x)` functions returns a `List` with the peptides
##'   that defined the connected components of object `x`.
##'
##' - The number of connected components of object `x` can be
##'   retrieved with `length(x)`.
##'
##' - The size of the connected components of object `x`, i.e the
##'   number of proteins in each component, can be retrieved with
##'   `lengths(x)`.
##'
##' - The `connectedComponents(x, i, simplify = TRUE)` function
##'   returns the peptide-by-protein sparse matrix (or `List` of
##'   matrices, if `length(i) > 1`), i.e. the subset of the adjacency
##'   matrix defined by the proteins in connected component(s)
##'   `i`. `i` is the numeric index (between 1 and `length(x)`) of the
##'   connected connected. If simplify is `TRUE` (default), then a
##'   matrix is returned instead of a `List` of matrices of length
##'   1. If set to `FALSE`, a `List` is always returned, irrespective
##'   of its length.
##'
##' @examples
##'
##' ## --------------------------------
##' ## From an adjacency matrix
##' ## --------------------------------
##' library(Matrix)
##' adj <- sparseMatrix(i = c(1, 2, 3, 3, 4, 4, 5),
##'                     j = c(1, 2, 3, 4, 3, 4, 5),
##'                     x = 1,
##'                     dimnames = list(paste0("Pep", 1:5),
##'                                    paste0("Prot", 1:5)))
##' adj
##' cc <- ConnectedComponents(adj)
##' cc
##'
##' length(cc)
##' lengths(cc)
##'
##' adjacencyMatrix(cc) ## same as adj
##' ccMatrix(cc)
##' ccList(cc)
##' ccPeptides(cc)
##' connectedComponents(cc, 3) ## a singel matrix
##' connectedComponents(cc, 1:2) ## a List
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
##'
##' length(cc)
##' table(lengths(cc))
##'
##' (i <- which(lengths(cc)  == 4))
##' ccomp <- connectedComponents(cc, i)
##'
##' ## A group of 4 proteins that all share peptide RTRYQAEVR
##' ccomp[[1]]
##'
##' ## A group of 4 proteins formed by 7 peptides: THPAERKPRRRKKR is
##' ## found in the two first proteins, KPTARRRKRK was found twice in
##' ## ECA3389, VVPVGLRALVWVQR was found in all 4 proteins, KLKPRRR
##' ## is specific to ECA3399, ...
##' ccomp[[2]]
NULL

setClass("ConnectedComponents",
         slots = c(adjMatrix = "Matrix",
                   ccMatrix = "Matrix",
                   ccList = "List",
                   ccPeptides = "List"))

##' @importFrom methods new
##'
##' @importFrom Matrix t tcrossprod rowSums
##'
##' @rdname ConnectedComponents
##'
##' @export
##'
##' @param object For the `ConnectedComponents` class constructor,
##'     either a sparse adjacency matrix of class `Matrix` or an
##'     instance of class `PSM`.
ConnectedComponents <- function(object) {
    if (is(object, "PSM")) {
        adj <- makeAdjacencyMatrix(object)
    } else if (is(object, "Matrix")) {
        adj <- object
    } else stop("'object' must be of class 'PSM' or 'Matrix.")
    if (is.null(colnames(adj)) | is.null(rownames(adj)))
        stop("The adjacency matrix used to create the object must have row and column names.")
    getCCpeptides <- function(cc, adj) {
        res <- adj[, cc, drop = FALSE]
        res <- res[Matrix::rowSums(res) > 0, , drop = FALSE]
        res
    }
    cc <- Matrix::tcrossprod(t(adj))
    n <- ncol(cc)
    cc_pep <- cc_list <- vector("list", length = n)
    i <- 1
    while (i <= n) {
        j <- i:n
        k <- which(cc[i, j] != 0)
        cc_list[[i]] <- colnames(adj)[j][k]
        cc_pep[[i]] <- rownames(getCCpeptides(cc_list[[i]], adj))
        i <- i + length(k)
    }
    sel <- lengths(cc_list) > 0
    new("ConnectedComponents",
        adjMatrix = adj,
        ccMatrix = cc,
        ccList = List(cc_list[sel]),
        ccPeptides = List(cc_pep[sel]))
}


setMethod("show", "ConnectedComponents",
          function(object) {
              cat(sprintf("An instance of class %s", class(object)), "\n")
              cat(" Number of proteins:", nrow(object@ccMatrix), "\n")
              cat(" Number of components:", length(object), "\n")
              cat(" Number of components by size:\n")
              tab <- table(lengths(object))
              msg <- strwrap(paste(paste0(tab, "(", names(tab), ")"),
                                   collapse = " "))
              message(paste(" ", msg, collapse = "\n"))
          })

##' @importFrom ProtGenerics adjacencyMatrix
##'
##' @export
##'
##' @rdname ConnectedComponents
setMethod("adjacencyMatrix", "ConnectedComponents",
          function(object) object@adjMatrix)

##' @export
##'
##' @rdname ConnectedComponents
##'
##' @param x An object of class `ConnectedComponents`.
ccMatrix <- function(x) {
    stopifnot(is(x, "ConnectedComponents"))
    x@ccMatrix
}

##' @export
##'
##' @rdname ConnectedComponents
ccPeptides <- function(x) {
    stopifnot(is(x, "ConnectedComponents"))
    x@ccPeptides
}

##' @export
##'
##' @rdname ConnectedComponents
ccList <- function(x) {
    stopifnot(is(x, "ConnectedComponents"))
    x@ccList
}

##' @export
##'
##' @rdname ConnectedComponents
connectedComponents <- function(x, i, simplify = TRUE) {
    if (any(i > length(x@ccList)))
        stop("Subscript out of bounds.")
    ans <- lapply(i, function(ii)
        x@adjMatrix[x@ccPeptides[[ii]],
                    x@ccList[[ii]],
                    drop = FALSE])
    if (length(ans) == 1 & simplify)
        return(ans[[1]])
    else return(List(ans))
}

##' ##' @export
##'
##' @rdname ConnectedComponents
setMethod("length", "ConnectedComponents",
          function(x) length(x@ccList))

##' ##' @export
##'
##' @rdname ConnectedComponents
setMethod("lengths", "ConnectedComponents",
          function(x) lengths(x@ccList))

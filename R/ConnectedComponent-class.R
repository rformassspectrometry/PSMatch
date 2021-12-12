##' @title Connected components
##'
##' @name ConnectedComponents
##'
##' @aliases ConnectedComponents ConnectedComponents-class ConnectedComponents length,ConnectedComponents lengths,ConnectedComponents adjacencyMatrix,ConnectedComponents ccMatrix show,ConnectedComponents connectedComponents [,ConnectedComponents,numeric,ANY,ANY [,ConnectedComponents,logical,ANY,ANY [,ConnectedComponents,integer,ANY,ANY
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
##' @slot adjMatrices A `List` containing adjacency matrices of each
##'     connected components.
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
##' - The number of connected components of object `x` can be
##'   retrieved with `length(x)`.
##'
##' - The size of the connected components of object `x`, i.e the
##'   number of proteins in each component, can be retrieved with
##'   `lengths(x)`.
##'
##' - The `connectedComponents(x, i, simplify = TRUE)` function
##'   returns the peptide-by-protein sparse adjacency matrix (or
##'   `List` of matrices, if `length(i) > 1`), i.e. the subset of the
##'   adjacency matrix defined by the proteins in connected
##'   component(s) `i`. `i` is the numeric index (between 1 and
##'   `length(x)`) of the connected connected. If simplify is `TRUE`
##'   (default), then a matrix is returned instead of a `List` of
##'   matrices of length 1. If set to `FALSE`, a `List` is always
##'   returned, irrespective of its length.
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
##' adjacencyMatrix(cc) ## same as adj above
##' ccMatrix(cc)
##'
##' connectedComponents(cc)
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
##' (i <- which(lengths(cc) == 4))
##' ccomp <- connectedComponents(cc, i)
##'
##' ## A group of 4 proteins that all share peptide RTRYQAEVR
##' ccomp[[1]]
##'
##' ## Visualise the adjacency matrix - here, we see how the single
##' ## peptides 'unites' the fous proteins.
##' plotAdjacencyMatrix(ccomp[[1]])
##'
##' ## A group of 4 proteins formed by 7 peptides: THPAERKPRRRKKR is
##' ## found in the two first proteins, KPTARRRKRK was found twice in
##' ## ECA3389, VVPVGLRALVWVQR was found in all 4 proteins, KLKPRRR
##' ## is specific to ECA3399, ...
##' ccomp[[3]]
##'
##' ## See how VVPVGLRALVWVQR is shared by ECA3406 ECA3415 ECA3389 and
##' ## links the three other componennts, namely ECA3399, ECA3389 and
##' ## (ECA3415, ECA3406). Filtering that peptide out would split that
##' ## protein group in three.
##' plotAdjacencyMatrix(ccomp[[3]])
##'
##' ## To select non-trivial components of size > 1
##' cc2 <- cc[lengths(cc) > 1]
##' cc2
NULL

setClass("ConnectedComponents",
         slots = c(adjMatrix = "Matrix",
                   ccMatrix = "Matrix",
                   adjMatrices = "List"))

##' @importFrom methods new
##'
##' @importFrom Matrix t tcrossprod rowSums
##'
##' @importFrom igraph graph_from_adjacency_matrix components groups
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
    cc <- Matrix::tcrossprod(t(adj))
    clu <- components(graph_from_adjacency_matrix(cc))
    adj_matrices <- lapply(unname(groups(clu)),
                           function(x) {
                               ans <- adj[, x, drop = FALSE]
                               ans[Matrix::rowSums(ans) > 0, , drop = FALSE]
                           })
    new("ConnectedComponents",
        adjMatrix = adj,
        ccMatrix = cc,
        adjMatrices = List(adj_matrices))
}


setMethod("show", "ConnectedComponents",
          function(object) {
              cat(sprintf("An instance of class %s", class(object)), "\n")
              cat(" Number of proteins:", nrow(object@ccMatrix), "\n")
              cat(" Number of components:", length(object@adjMatrices), "\n")
              cat(" Number of components by size:\n")
              tab <- table(lengths(object))
              msg <- strwrap(paste(paste0(tab, "(", names(tab), ")"),
                                   collapse = " "))
              message(paste(" ", msg, collapse = "\n"))
          })


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
##'
##' @param i `integer()` with the index of the component(s) to return.
##'
##' @param simplify `logical(1)` if `TRUE` (default), the output is
##'     simplified to sparse matrix if `i` was of length 1, otherwise
##'     a `List` is returned. Always a `List` if `FALSE`.
connectedComponents <- function(x, i, simplify = TRUE) {
    if (missing(i))
        return(x@adjMatrices)
    stopifnot(is.numeric(i))
    if (any(i > length(x)))
        stop("Subscript out of bounds.")
    ans <- x@adjMatrices[i]
    if (length(ans) == 1 & simplify)
        return(ans[[1]])
    else return(ans)
}

plotAdjacencyMatrix <- function(x, layout = layout_nicely) {
    g <- graph_from_incidence_matrix(x)
    V(g)$color <- ifelse(names(V(g)) %in% colnames(x), "steelblue", "orange")
    plot(g, layout = layout_nicely)
}

##' @export
##'
##' @rdname ConnectedComponents
setMethod("length", "ConnectedComponents",
          function(x) length(x@adjMatrices))

##' ##' @export
##'
##' @rdname ConnectedComponents
setMethod("lengths", "ConnectedComponents",
          function(x) sapply(x@adjMatrices, ncol))

##' @export
##'
##' @rdname ConnectedComponents
##'
##' @param i `numeric()`, `integer()` or `logical()` to subset the
##'     `ConnectedComponents` instance. If a `logical()`, it must be
##'     of same length as the object is subsets.
##'
##' @param j ignored
##'
##' @param ... ignored
##'
##' @param drop ignore
setMethod("[", c("ConnectedComponents", "integer"),
          function(x, i, j, ..., drop = FALSE) {
    if (!missing(j))
        stop("Subsetting ConnectedComponents by 'i' only.")
    if (missing(i))
        return(x)
    subsetConnectedComponents(x, as.integer(i))
})

##' @export
##'
##' @rdname ConnectedComponents
setMethod("[", c("ConnectedComponents", "logical"),
          function(x, i, j, ..., drop = FALSE) {
              if (length(i) != length(x))
                  stop("'i' must be of same length than 'x'.")
              x[which(i)]
          })

##' @export
##'
##' @rdname ConnectedComponents
setMethod("[", c("ConnectedComponents", "numeric"),
          function(x, i, j, ..., drop = FALSE) x[as.integer(i)])

##' @importFrom Matrix colSums
subsetConnectedComponents <- function(object, i) {
    stopifnot(is.integer(i))
    stopifnot(is(object, "ConnectedComponents"))
    stopifnot(max(i) <= length(object))
    ## peptides in connected components to be kept
    keep_peptides <- unique(unlist(lapply(object@adjMatrices[i], rownames)))
    adj <- object@adjMatrix[keep_peptides, , drop = FALSE]
    ## remove protein without any peptides
    adj <- adj[, Matrix::colSums(adj) > 0, drop = FALSE]
    ## New object from updated adjacency matrix
    ConnectedComponents(adj)
}

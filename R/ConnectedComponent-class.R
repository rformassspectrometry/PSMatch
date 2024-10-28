##' @title Connected components
##'
##' @name ConnectedComponents
##'
##' @aliases ConnectedComponents ConnectedComponents-class ConnectedComponents length,ConnectedComponents adjacencyMatrix,ConnectedComponents ccMatrix show,ConnectedComponents connectedComponents [,ConnectedComponents,numeric,ANY,ANY [,ConnectedComponents,logical,ANY,ANY [,ConnectedComponents,integer,ANY,ANY dims,ConnectedComponents ncols,ConnectedComponents nrows,ConnectedComponents prioritiseConnectedComponents prioritizeConnectedComponents
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
##'   `ncols(x)`. The number of peptides defining the connected
##'   components can be retrieved with `nrows(x)`. Both can be
##'   accessed with `dims(x)`.
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
##' - To help with the exploration of individual connected Components,
##'   the `prioritiseConnectedComponents()` function will take an
##'   instance of `ConnectedComponents` and return a `data.frame` where
##'   the component indices are ordered based on their potential to
##'   clean up/flag some peptides and split protein groups in small
##'   groups or individual proteins, or simply explore them. The
##'   prioritisation is based on a set of metrics computed from the
##'   component's adjacency matrix, including its dimensions, row and
##'   col sums maxima and minima, its sparsity and the number of
##'   communities and their modularity that quantifies how well the
##'   communities separate (see [modularity.igraph()]. Note that
##'   trivial components, i.e. those composed of a single peptide and
##'   protein are excluded from the prioritised results. This
##'   `data.frame` is ideally suited for a principal component
##'   analysis (using for instance [prcomp()]) for further inspection
##'   for component visualisation with [plotAdjacencyMatrix()].
##'
##' @return The `ConnectedComponents()` constructor returns an
##'     instance of class `ConnectedComponents`. The *Creating and
##'     manipulating objects* section describes the return values of
##'     the functions that manipulate `ConnectedComponents` objects.
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
##' ncols(cc)
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
##' table(ncols(cc))
##'
##' (i <- which(ncols(cc) == 4))
##' ccomp <- connectedComponents(cc, i)
##'
##' ## A group of 4 proteins that all share peptide RTRYQAEVR
##' ccomp[[1]]
##'
##' ## Visualise the adjacency matrix - here, we see how the single
##' ## peptides (white node) 'unites' the four proteins (blue nodes)
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
##' ## Colour protein node based on protein names similarity
##' plotAdjacencyMatrix(ccomp[[3]], 1)
##'
##' ## To select non-trivial components of size > 1
##' cc2 <- cc[ncols(cc) > 1]
##' cc2
##'
##' ## Use components features to prioritise their exploration
##' pri_cc <- prioritiseConnectedComponents(cc)
##' pri_cc
##'
##' plotAdjacencyMatrix(connectedComponents(cc, 1082), 1)
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
##'
##' @param ... Additional arguments passed to
##'     [`makeAdjacencyMatrix()`] when `object` is of class [PSM()].
ConnectedComponents <- function(object, ...) {
    if (is(object, "PSM")) {
        adj <- makeAdjacencyMatrix(object, ...)
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
    ## Set the order of the CC adjacency matrix list using the name of
    ## the first protein, otherwise it will be randomly defined by
    ## igraph::components/groups, which is very annoying when
    ## analysing data.
    o <- order(vapply(adj_matrices, function(x) colnames(x)[1], ""))
    new("ConnectedComponents",
        adjMatrix = adj,
        ccMatrix = cc,
        adjMatrices = List(adj_matrices[o]))
}


setMethod("show", "ConnectedComponents",
          function(object) {
              cat(sprintf("An instance of class %s", class(object)), "\n")
              cat(" Number of proteins:", nrow(ccMatrix(object)), "\n")
              cat(" Number of components:", length(connectedComponents(object)), "\n")
              cat(" Number of components [peptide x proteins]:\n  ")
              dim_mat <- dims(object)
              cat(paste0(sum(dim_mat[, 1] == 1 & dim_mat[, 2] == 1), "[1 x 1] "))
              cat(paste0(sum(dim_mat[, 1] == 1 & dim_mat[, 2] > 1), "[1 x n] "))
              cat(paste0(sum(dim_mat[, 1] > 1 & dim_mat[, 2] == 1), "[n x 1] "))
              cat(paste0(sum(dim_mat[, 1] > 1 & dim_mat[, 2] > 1), "[n x n]\n"))
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

##' @export
##'
##' @rdname ConnectedComponents
setMethod("length", "ConnectedComponents",
          function(x) length(connectedComponents(x)))

##' @export
##'
##' @importFrom BiocGenerics dims
##'
##' @rdname ConnectedComponents
setMethod("dims", "ConnectedComponents",
          function(x) {
              ans <- t(vapply(connectedComponents(x), dim, c(1, 1)))
              colnames(ans) <- c("nrow", "ncol")
              ans
          })

##' @export
##'
##' @importFrom BiocGenerics ncols
##'
##' @rdname ConnectedComponents
setMethod("ncols", "ConnectedComponents",
          function(x) vapply(connectedComponents(x), ncol, 1L))

##' @export
##'
##' @importFrom BiocGenerics nrows
##'
##' @rdname ConnectedComponents
setMethod("nrows", "ConnectedComponents",
          function(x) vapply(connectedComponents(x), nrow, 1L))

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
              if (!missing(j))
                  stop("Subsetting ConnectedComponents by 'i' only.")
              if (length(i) != length(x))
                  stop("'i' must be of same length than 'x'.")
              x[which(i)]
          })

##' @export
##'
##' @rdname ConnectedComponents
setMethod("[", c("ConnectedComponents", "numeric"),
          function(x, i, j, ..., drop = FALSE) {
              if (!missing(j))
                  stop("Subsetting ConnectedComponents by 'i' only.")
              x[as.integer(i)]
          })

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

##' @export
##'
##' @rdname ConnectedComponents
##'
##' @importFrom igraph cluster_louvain modularity
prioritiseConnectedComponents <- function(x) {
    stopifnot(is(x, "ConnectedComponents"))
    ans <- data.frame(dims(x))
    cc_x <- connectedComponents(x)
    ## community metrics
    com_metrics <- t(vapply(cc_x,
                     function(xx) {
                         g <- graph_from_biadjacency_matrix(xx)
                         com <- cluster_louvain(g)
                         c(n_coms = length(com),
                           mod_coms = modularity(com))
                     }, c(1, 1)))
    ans <- cbind(ans, com_metrics)
    ## adjacency matrix metrics
    adj_metrics <- t(vapply(cc_x,
                            function(xx)
                                c(n = sum(xx),
                                  rs_min = min(rowSums(xx)),
                                  rs_max = max(rowSums(xx)),
                                  cs_min = min(colSums(xx)),
                                  cs_max = max(colSums(xx)),
                                  sparsity = sum(xx == 0)/prod(dim(xx))),
                            numeric(6)))
    ans <- cbind(ans, adj_metrics)
    sel <- ans$ncol > 1 & ans$nrow > 1
    ans <- ans[sel, ]
    ans[order(ans$mod_coms, decreasing = TRUE), ]
}

##' @export
##'
##' @rdname ConnectedComponents
prioritizeConnectedComponents <- prioritiseConnectedComponents

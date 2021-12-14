##' @title Describe protein and peptide compositions
##'
##' @description
##'
##' It is important to explore PSM results prior to any further
##' downstream analysies. Two functions, that work on [PSM()] and
##' [ConnectedComponents()] objects can be used for this:
##'
##' - The `describeProteins()` function describe protein composition
##'   in terms of unique and shared peptides.
##'
##' - The `describePeptides()` function describe unique/shared peptide
##'   composition.
##'
##' @name describeProteins
##'
##' @aliases describePeptides
##'
##' @rdname describeProteins
##'
##' @examples
##' f <- msdata::ident(full.names = TRUE, pattern = "TMT")
##' basename(f)
##' psm <- PSM(f) |>
##'        filterPsmDecoy() |>
##'        filterPsmRank()
##'
##' describePeptides(psm)
##' describeProteins(psm)
NULL

## --------------------------------------------------------------------
## Describe protein composition in terms of unique and shared peptides
## --------------------------------------------------------------------

##' @importFrom Matrix rowSums
.describeProteins <- function(adj) {
    unique_peps <- names(which(Matrix::rowSums(adj) == 1))
    ans_unique <- ans_shared <- rep(NA, ncol(adj))
    for (i in seq_along(ans_unique)) {
        peps_i <- names(which(adj[, i] != 0))
        ans_unique[i] <- all(peps_i %in% unique_peps)
        ans_shared[i] <- all(!peps_i %in% unique_peps)
    }
    ans_mixed <- !(ans_unique | ans_shared)
    ans <- data.frame(uniqueOnly = ans_unique,
                      sharedOnly = ans_shared,
                      both = ans_mixed,
                      row.names = colnames(adj))
    message(ncol(adj), " proteins composed of")
    message(" only unique peptides: ", sum(ans$uniqueOnly))
    message(" only shared peptides: ", sum(ans$sharedOnly))
    message(" unique and shared peptides: ", sum(ans$both))
    invisible(ans)
}


##' @export
##'
##' @importFrom methods is
##'
##' @rdname describeProteins
##'
##' @param object Either an instance of class `Matrix`, [PSM()] or
##'     [ConnectedComponents()].
describeProteins <- function(object) {
    if (is(object, "PSM")) {
        adj <- makeAdjacencyMatrix(object)
    } else if (is(object, "ConnectedComponents")) {
        adj <- adjacencyMatrix(object)
    } else if (is(object, "Matrix")) {
        adj <- object
    } else stop("Object must be of class 'Matrix', 'PSM' or 'ConnectedComponents'")
    .describeProteins(adj)
}

## -------------------------------------------
## Describe unique/shared peptide composition
## -------------------------------------------
.describePeptides <- function(adj) {
    tab <- table(Matrix::rowSums(adj))
    message(nrow(adj), " peptides composed of")
    message(" unique peptides: ", tab["1"])
    message(" shared peptides (among protein):")
    tab2 <- tab[names(tab) != "1"]
    msg <- strwrap(paste(paste0(tab2, "(", names(tab2), ")"),
                         collapse = " "))
    message(paste(" ", msg, collapse = "\n"))
    invisible(tab)
}

##' @export
##'
##' @rdname describeProteins
describePeptides <- function(object) {
    if (is(object, "PSM")) {
        adj <- makeAdjacencyMatrix(object)
    } else if (is(object, "ConnectedComponents")) {
        adj <- adjacencyMatrix(object)
    } else if (is(object, "Matrix")) {
        adj <- object
    } else stop("Object must be of class 'Matrix', 'PSM' or 'ConnectedComponents'")
 .describePeptides(adj)
}

## --------------------------------------------------------------------
## Describe protein composition in terms of unique and shared peptides
## --------------------------------------------------------------------

.describeProteins <- function(adj) {
    unique_peps <- names(which(rowSums(adj) == 1))
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

describeProteins <- function(object) {
    if (is(object, "PSM")) {
        if (is.null(metadata(object)$adjacencyMatrix))
            metadata(object)$adjacencyMatrix <- makeAdjacencyMatrix(object)
        ans <- .describeProteins(metadata(object)$adjacencyMatrix)
    } else if (is(object, "ConnectedComponents")) {
        ans <- .describeProteins(object@adjacencyMatrix)
    } else stop("Object must be of class 'PSM' or 'ConnectedComponents'")
    ans
}

## -------------------------------------------
## Describe unique/shared peptide composition
## -------------------------------------------

.describePeptides <- function(adj) {
    tab <- table(rowSums(adj))
    message(nrow(adj), " peptides composed of")
    message(" unique peptides: ", tab[1])
    message(" shared peptides (among protein):")
    msg <- strwrap(paste(paste0(tab[-1], "(", names(tab[-1]), ")"),
                         collapse = " "))
    message(paste(" ", msg, collapse = "\n"))
    invisible(tab)
}

describePeptides <- function(object) {
    if (is(object, "PSM")) {
        if (is.null(metadata(object)$adjacencyMatrix))
            metadata(object)$adjacencyMatrix <- makeAdjacencyMatrix(object)
        ans <- .describePeptides(metadata(object)$adjacencyMatrix)
    } else if (is(object, "ConnectedComponents")) {
        ans <- .describePeptides(object@adjacencyMatrix)
    } else stop("Object must be of class 'PSM' or 'ConnectedComponents'")
    ans
}

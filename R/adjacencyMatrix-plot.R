##' @importFrom igraph layout_nicely graph_from_incidence_matrix V V<- plot.igraph
##'
##' @export
##'
##' @importFrom utils adist
##'
##' @name adjacencyMatrix
##'
##' @param x A peptide-by-protein adjacency matrix.
##'
##' @param colLevel `numeric(1)` indicating the protein colouring
##'     level to use. If 0 (default), all protein nodes are labelled
##'     in steelblue. For values > 0, approximate string distances
##'     (see [adist()] between protein names are calculated and nodes
##'     of protein that have names that differ will be coloured
##'     differently, with higher values leading to more colours. While
##'     no maximum to this value is defined in the code, it shouldn't
##'     be higher than the number of proteins. Peptide nodes are
##'     always labelled in orange.
##'
##' @param layout A graph layout, as defined in the `ipgraph`
##'     package. Default is [layout_as_bipartite()].
plotAdjacencyMatrix <- function(x, colLevel = 0,
                                layout = igraph::layout_nicely) {
    g <- graph_from_incidence_matrix(x)
    colLevel <- colLevel[1]
    if (colLevel == 0) {
        V(g)$color <- ifelse(names(V(g)) %in% colnames(x), "steelblue", "orange")
    } else {
        ## distances between protein names
        pnames0 <- colnames(x)
        ## parse swissprot identifiers to keep gene names, but keep
        ## _SPECIES suffix to take it into account proteins from
        ## different species. Identical suffixes don't affect the
        ## distances, so can be safely kept for withing species
        ## distances.
        pnames <- sub("^.+\\|", "", pnames0)
        pnm <- adist(pnames)
        rownames(pnm) <- colnames(pnm) <- pnames
        for (i in seq_len(colLevel))
            pnm[pnm == max(pnm)] <- 0
        ## set colours from clique memberships
        cls <- components(graph_from_adjacency_matrix(pnm))$membership + 1
        V(g)$color <- rep(1, length(V(g)))
        V(g)$color[match(pnames0, names(V(g)))] <- cls
    }
    plot.igraph(g, layout = layout)
    invisible(g)
}

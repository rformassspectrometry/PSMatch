##' @importFrom igraph layout_nicely graph_from_incidence_matrix V V<- plot.igraph
##'
##' @export
##'
##' @name adjacencyMatrix
##'
##' @param x A peptide-by-protein adjacency matrix.
##'
##' @param layout A graph layout, as defined in the `ipgraph`
##'     package. Default is [layout_as_bipartite()].
plotAdjacencyMatrix <- function(x, layout = igraph::layout_as_bipartite) {
    g <- graph_from_incidence_matrix(x)
    V(g)$color <- ifelse(names(V(g)) %in% colnames(x), "steelblue", "orange")
    plot.igraph(g, layout = layout_nicely)
    invisible(g)
}

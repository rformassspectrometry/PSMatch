##' @importFrom igraph layout_nicely graph_from_biadjacency_matrix V V<- plot.igraph set_vertex_attr
##'
##' @export
##'
##' @importFrom utils adist
##'
##' @name adjacencyMatrix
##'
##' @param protColors Either a `numeric(1)` or a named `character()`
##'     of colour names. The numeric value indicates the protein
##'     colouring level to use. If 0 (default), all protein nodes are
##'     labelled in steelblue. For values > 0, approximate string
##'     distances (see [adist()]) between protein names are calculated
##'     and nodes of proteins that have names that differ will be
##'     coloured differently, with higher values leading to more
##'     colours. While no maximum to this value is defined in the
##'     code, it shouldn't be higher than the number of proteins. If a
##'     character is used, it should be a character of colour names
##'     named by protein identifiers. That vector should provide
##'     colours for at least all proteins in the adjacency matrix `m`,
##'     but more protein could be named. The latter is useful when
##'     generating a colour vector for all proteins in a dataset and
##'     use it for different adjacency matrix visualisations.
##'
##' @param pepColors Either `NULL` (default) for no peptide colouring
##'     (white nodes) or a named `character()` of colour names. It
##'     should be a character of colour names named by peptide
##'     identifiers. That vector should provide colours for at least
##'     all peptides in the adjacency matrix `m`, but more peptides
##'     could be named. The latter is useful when generating a colour
##'     vector for all peptides in a dataset and use it for different
##'     adjacency matrix visualisations.
##'
##' @param layout A graph layout, as defined in the `ipgraph`
##'     package. Default is [igraph::layout_as_bipartite()].
plotAdjacencyMatrix <- function(m,
                                protColors = 0,
                                pepColors = NULL,
                                layout = igraph::layout_nicely) {
    g <- graph_from_biadjacency_matrix(m)
    if (is.character(protColors)) {
        ## Expecting a named vector of colour characters
        if (is.null(names(protColors)))
            stop("Protein colours must be named.")
        ## These names must match the proteins/columns in the
        ## input adjacency matrix.
        if (!all(colnames(m) %in% names(protColors)))
            stop("Please provide colours for all the proteins in the adjacency matrix.")
        ## Keeping relevant protein colours - there could be more, if
        ## that vector was generated for all proteins in the dataset,
        ## and not specifically for this adjacency matrix.
        protColors <- protColors[colnames(m)]
        V(g)$color[match(names(protColors), names(V(g)))] <- protColors
    } else if (is.numeric(protColors)) {
        ## Setting colours automatically based - single colour for all
        ## proteins in steelblue
        proColors <- protColors[1]
        if (proColors == 0) {
            V(g)$color <- ifelse(names(V(g)) %in% colnames(m),
                                 "steelblue", NA_character_)
        } else {
            ## Set colour based on the distances between protein names
            pnames0 <- colnames(m)
            ## Parse swissprot identifiers to keep gene names, but
            ## keep _SPECIES suffix to take into account proteins from
            ## different species. Identical suffixes don't affect the
            ## distances, so can be safely kept for same-species
            ## distances.
            pnames <- sub("^.+\\|", "", pnames0)
            pnm <- adist(pnames)
            rownames(pnm) <- colnames(pnm) <- pnames
            for (i in seq_len(protColors))
                pnm[pnm == max(pnm)] <- 0
            ## set colours from clique memberships
            cls <- components(graph_from_adjacency_matrix(pnm))$membership + 1
            V(g)$color <- rep(0, length(V(g)))
            V(g)$color[match(pnames0, names(V(g)))] <- cls
        }
    } else stop("protColors must be a named character or a single numeric")
    if (!is.null(pepColors)) {
        ## Expecting a named vector of colour characters
        if (!is.character(pepColors) || is.null(names(pepColors)))
            stop("pepColors must be a named character or NULL.")
        ## These names must match the peptides/rows in the input
        ## adjacency matrix.
        if (!all(rownames(m) %in% names(pepColors)))
            stop("Please provide colours for all the peptides in the adjacency matrix.")
        ## Keeping relevant protein colours - there could be more, if
        ## that vector was generated for all proteins in the dataset,
        ## and not specifically for this adjacency matrix.
        pepColors <- pepColors[rownames(m)]
        V(g)$color[match(names(pepColors), names(V(g)))] <- pepColors
    }
    pep_nodes <- which(names(V(g)) %in% rownames(m))
    prot_nodes <- which(names(V(g)) %in% colnames(m))
    v_shapes <- rep("none", sum(dim(m)))
    v_shapes[pep_nodes] <- "circle"
    v_shapes[prot_nodes] <- "square"
    g <- igraph::set_vertex_attr(g, "shape", value = v_shapes)
    plot.igraph(g, layout = layout)
    invisible(g)
}
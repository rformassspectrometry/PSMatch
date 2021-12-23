test_that("ConnectedComponents works", {
    ans <- ConnectedComponents(adj)
    expect_null(show(ans))
    expect_identical(adj, adjacencyMatrix(ans))
    expect_equivalent(cc, ccMatrix(ans))
    expect_identical(dimnames(cc), dimnames(ccMatrix(ans)))
    expect_identical(adjMatrices[[1]],
                     connectedComponents(ans, 1))
    expect_identical(length(ans), 4L)
    expect_identical(length(connectedComponents(ans)), 4L)
    expect_identical(length(connectedComponents(ans, 1:3)), 3L)
    expect_identical(length(connectedComponents(ans, 1, simplify = FALSE)), 1L)
    expect_identical(connectedComponents(ans, 1), adjMatrices[[1]])
    expect_identical(ncols(ans), c(1L, 1L, 2L, 3L))
    expect_identical(nrows(ans), c(2L, 1L, 2L, 3L))
    expect_identical(connectedComponents(ans[1:2]), adjMatrices[1:2])
    expect_error(connectedComponents(ans, 5), "Subscript out of bounds.")
    expect_error(ans[TRUE], "'i' must be of same length than 'x'.")
    expect_identical(ans[c(TRUE, FALSE, FALSE, FALSE)], ans[1])
    expect_error(ans[i = 1, j = 1],
                 "Subsetting ConnectedComponents by 'i' only.")
    expect_error(ans[i = 1L, j = 1],
                 "Subsetting ConnectedComponents by 'i' only.")
    expect_error(ans[c(TRUE, FALSE, FALSE, FALSE), j = 1],
                 "Subsetting ConnectedComponents by 'i' only.")
})


test_that("ConnectedComponents works from PSM", {
    cc1 <- ConnectedComponents(adj)
    cc2 <- ConnectedComponents(psmdf)
    cc2@adjMatrix <- cc2@adjMatrix[order(rownames(cc2@adjMatrix)), ]
    expect_identical(cc1, cc2)
})



test_that("prioritiseConnectedComponents() works", {
    cc <- ConnectedComponents(adj)
    p1 <- prioritiseConnectedComponents(cc)
    p2 <- prioritizeConnectedComponents(cc)
    expect_identical(p1, p2)
    ## Check CC 4
    cc_i <- connectedComponents(cc, 4)
    expect_identical(p1["4", "ncol"], ncol(cc_i))
    expect_identical(p1["4", "nrow"], nrow(cc_i))
    expect_identical(p1["4", "n"], sum(cc_i))
    expect_identical(p1["4", "rs_min"], min(rowSums(cc_i)))
    expect_identical(p1["4", "rs_max"], max(rowSums(cc_i)))
    expect_identical(p1["4", "cs_min"], min(colSums(cc_i)))
    expect_identical(p1["4", "cs_max"], max(colSums(cc_i)))
    expect_identical(p1["4", "sparsity"], sum(cc_i == 0)/(ncol(cc_i) * nrow(cc_i)))
    cl <- cluster_louvain(igraph::graph_from_incidence_matrix(cc_i))
    expect_identical(p1["4", "n_coms"], as.numeric(length(cl)))
    expect_identical(p1["4", "mod_coms"], modularity(cl))
    ## Check CC 3
    cc_i <- connectedComponents(cc, 3)
    expect_identical(p1["3", "ncol"], ncol(cc_i))
    expect_identical(p1["3", "nrow"], nrow(cc_i))
    expect_identical(p1["3", "n"], sum(cc_i))
    expect_identical(p1["3", "rs_min"], min(rowSums(cc_i)))
    expect_identical(p1["3", "rs_max"], max(rowSums(cc_i)))
    expect_identical(p1["3", "cs_min"], min(colSums(cc_i)))
    expect_identical(p1["3", "cs_max"], max(colSums(cc_i)))
    expect_identical(p1["3", "sparsity"], sum(cc_i == 0)/(ncol(cc_i) * nrow(cc_i)))
    cl <- cluster_louvain(igraph::graph_from_incidence_matrix(cc_i))
    expect_identical(p1["3", "n_coms"], as.numeric(length(cl)))
    expect_identical(p1["3", "mod_coms"], modularity(cl))
})

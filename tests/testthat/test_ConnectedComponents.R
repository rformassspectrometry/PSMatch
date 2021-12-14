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
    expect_identical(lengths(ans), c(1L, 1L, 2L, 3L))
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

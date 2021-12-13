adj <- sparseMatrix(
    i = c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8),
    j = c(1, 2, 3, 4, 3, 4, 5, 6, 6, 7, 7, 1),
    x = 1,
    dimnames = list(
        paste0("p", 1:8),
        paste0("P", 1:7)))

cc <- sparseMatrix(
    i = c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7),
    j = c(1, 2, 3, 4, 3, 4, 5, 6, 5, 6, 7, 6, 7),
    x = c(2, 1, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 2),
    dimnames = list(
        paste0("P", 1:7),
        paste0("P", 1:7)))

adjMatrices <-
    List(
        sparseMatrix(i = c(1, 2), j = c(1, 1), x = 1,
                     dimnames = list(c("p1", "p8"),
                                     "P1")),
        sparseMatrix(i = 1, j = 1, x = 1,
                     dimnames = list("p2", "P2")),
        sparseMatrix(i = c(1, 1, 2, 2),
                     j = c(1, 2, 1, 2),
                     x = 1,
                     dimnames = list(
                         c("p3", "p4"),
                         c("P3", "P4"))),
        sparseMatrix(i = c(1, 1, 2, 2, 3),
                     j = c(1, 2, 2, 3, 3),
                     x = 1,
                     dimnames = list(
                         c("p5", "p6", "p7"),
                         c("P5", "P6", "P7")))
    )



test_that("ConnectedComponents works", {
    ans <- ConnectedComponents(adj)
    expect_identical(adj, adjacencyMatrix(ans))
    expect_equivalent(cc, ccMatrix(ans))
    expect_identical(dimnames(cc), dimnames(ccMatrix(ans)))
    expect_identical(adjMatrices[[1]],
                     connectedComponents(ans, 1))
})

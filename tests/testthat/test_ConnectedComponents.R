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
                         c("P5", "P6", "P7"))))


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




df <- data.frame(peptide = paste0("p", c(1, 8, 2, 3, 4, 3, 4, 5, 5, 6, 6, 7)),
           protein = rep(colnames(adj), colSums(adj)))
psmdf <- PSM(df, protein = "protein", peptide = "peptide")

test_that("ConnectedComponents works from PSM", {
    cc1 <- ConnectedComponents(adj)
    cc2 <- ConnectedComponents(psmdf)
    cc2@adjMatrix <- cc2@adjMatrix[order(rownames(cc2@adjMatrix)), ]
    expect_identical(cc1, cc2)
})

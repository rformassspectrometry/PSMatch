# setting R_TESTS to empty string because of
# https://github.com/hadley/testthat/issues/144
# revert this when that issue in R is fixed.
Sys.setenv("R_TESTS" = "")

library("testthat")
library("PSMatch")

adj <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8),
    j = c(1, 2, 3, 4, 3, 4, 5, 6, 6, 7, 7, 1),
    x = 1,
    dimnames = list(
        paste0("p", 1:8),
        paste0("P", 1:7)))

cc <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7),
    j = c(1, 2, 3, 4, 3, 4, 5, 6, 5, 6, 7, 6, 7),
    x = c(2, 1, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 2),
    dimnames = list(
        paste0("P", 1:7),
        paste0("P", 1:7)))

adjMatrices <-
    S4Vectors::List(
        Matrix::sparseMatrix(
                    i = c(1, 2), j = c(1, 1), x = 1,
                    dimnames = list(c("p1", "p8"),
                                    "P1")),
        Matrix::sparseMatrix(
                    i = 1, j = 1, x = 1,
                    dimnames = list("p2", "P2")),
        Matrix::sparseMatrix(
                    i = c(1, 1, 2, 2),
                    j = c(1, 2, 1, 2),
                    x = 1,
                    dimnames = list(
                        c("p3", "p4"),
                        c("P3", "P4"))),
        Matrix::sparseMatrix(
                    i = c(1, 1, 2, 2, 3),
                    j = c(1, 2, 2, 3, 3),
                    x = 1,
                    dimnames = list(
                        c("p5", "p6", "p7"),
                        c("P5", "P6", "P7"))))

psmdf <- PSM(data.frame(peptide = paste0("p", c(1, 8, 2, 3, 4, 3, 4, 5, 5, 6, 6, 7)),
                        protein = rep(colnames(adj), Matrix::colSums(adj))),
             protein = "protein", peptide = "peptide")

test_check("PSMatch")

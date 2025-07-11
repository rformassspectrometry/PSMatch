test_that("describePeptides works", {
    d1 <- describePeptides(adj)
    d2 <- describePeptides(psmdf)
    expect_identical(d1, d2)
    ## > rowSums(adj)
    ## p1 p2 p3 p4 p5 p6 p7 p8
    ##  1  1  2  2  2  2  1  1
    expect_identical(d1, table(c(1, 1, 2, 2, 2, 2, 1, 1)))
})


test_that("describePeptides works with non-binary adj matrix", {
    ccx <- new("dgCMatrix",
               i = c(0L, 2L, 4L, 5L, 0L, 2L, 4L, 5L, 1L, 3L,
                     5L, 5L, 6L),
               p = c(0L, 4L, 8L, 11L, 13L),
               Dim = c(7L, 4L),
               Dimnames = list(
                   c("THPAERKPRRRKKR", "KPTARRRKRK",
                     "PLAQGGQLNRLSAIRGLFR", "RRKRKPDSLKK", "KPRRRK",
                     "VVPVGLRALVWVQR", "KLKPRRR"),
                   c("ECA3406", "ECA3415", "ECA3389", "ECA3399")),
               x = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1),
               factors = list())
    expect_identical(describePeptides(ccx),
                     table(c(1, 1, 1, 2, 2, 2, 4)))
})


test_that("describeProteins works", {
    d1 <- describeProteins(adj)
    d2 <- describeProteins(psmdf)
    expect_identical(d1, d2)
    ans <- data.frame(P1 = c(TRUE, FALSE, FALSE),
               P2 = c(TRUE, FALSE, FALSE),
               P3 = c(FALSE, TRUE, FALSE),
               P4 = c(FALSE, TRUE, FALSE),
               P5 = c(FALSE, TRUE, FALSE),
               P6 = c(FALSE, TRUE, FALSE),
               P7 = c(FALSE, FALSE, TRUE),
               row.names = c("uniqueOnly", "sharedOnly",
                             "both"))
    ans <- data.frame(t(ans))
    expect_identical(d1, ans)
})


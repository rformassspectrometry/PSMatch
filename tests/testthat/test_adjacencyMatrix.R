m <- matrix(0, nrow = 5, ncol = 4)
colnames(m) <- LETTERS[1:4]
rownames(m) <- letters[1:5]
m[1,1] <- m[1, 3] <-
    m[2, 2] <-
    m[3, 2] <-
    m[4, 3] <- m[4, 4] <-
    m[5, 1] <- m[5, 3] <- 1
m <- Matrix::Matrix(m)

vec <- c("A;C", "B", "B", "C;D", "A;C")
names(vec) <- letters[1:5]

test_that("Compare makeAjacendyMatrix() makePeptideProteinVector() outputs on test data", {
    m2 <- makeAdjacencyMatrix(vec)
    ## first check that row/colnames as same
    expect_identical(sort(colnames(m)),
                     sort(colnames(m2)))
    expect_identical(sort(rownames(m)),
                     sort(rownames(m2)))
    ## reorder col/rows before checking identify
    m2 <- m2[, colnames(m)]
    m2 <- m2[rownames(m), ]
    expect_identical(m, m2)
    ## expecting same order here
    vec2 <- makePeptideProteinVector(m)
    expect_identical(vec, vec2)
})

test_that("Check makeAjacendyMatrix() works without split", {
    vec <- LETTERS[1:3]
    names(vec) <- letters[1:3]
    adj0 <- Matrix::Matrix(diag(3))
    colnames(adj0) <- LETTERS[1:3]
    rownames(adj0) <- letters[1:3]
    adj1 <- makeAdjacencyMatrix(vec)
    adj2 <- makeAdjacencyMatrix(vec, split = NULL)
    expect_equivalent(adj0, adj2)
    expect_equivalent(adj1, adj2)
})

test_that("Check makeAjacendyMatrix() fails with wrong input", {
    expect_error(makeAdjacencyMatrix(data.frame()))
    expect_error(makeAdjacencyMatrix(list()))
})

test_that("makeAjacendyMatrix() works on PMS data (1)", {
    f <- msdata::ident(full.names = TRUE, pattern = "TMT")
    psm <- filterPSMs(PSM(f))
    adj <- makeAdjacencyMatrix(psm)
    n_pep <- length(unique(psm[[psmVariables(psm)["peptide"]]]))
    expect_identical(nrow(adj), n_pep)
    n_prot <- length(unique(psm[[psmVariables(psm)["protein"]]]))
    expect_identical(ncol(adj), n_prot)
    adj2 <- makeAdjacencyMatrix(psm, binary = TRUE)
    expect_true(all(as.vector(adj2) %in% c(0, 1)))
    expect_identical(dimnames(adj), dimnames(adj2))
})

test_that("makeAjacendyMatrix() works on PMS data (2)", {
    f <- msdata::ident(full.names = TRUE, pattern = "TMT")
    psm <- PSM(f)
    adj <- makeAdjacencyMatrix(psm)
    n_pep <- length(unique(psm[[psmVariables(psm)["peptide"]]]))
    expect_identical(nrow(adj), n_pep)
    n_prot <- length(unique(psm[[psmVariables(psm)["protein"]]]))
    expect_identical(ncol(adj), n_prot)
    adj2 <- makeAdjacencyMatrix(psm, binary = TRUE)
    expect_true(all(as.vector(adj2) %in% c(0, 1)))
    expect_identical(dimnames(adj), dimnames(adj2))
})

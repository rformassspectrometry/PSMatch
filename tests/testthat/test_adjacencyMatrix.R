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
    ## try again without psmVariables
    metadata(psm)$variables["protein"] <- "NOTVALID"
    expect_error(makeAdjacencyMatrix(psm))
    expect_error(adjacencyMatrix(psm))
    metadata(psm)$variables["protein"] <- "DatabaseAccess"
    metadata(psm)$variables["peptide"] <- "NOTVALID"
    expect_error(makeAdjacencyMatrix(psm))
    metadata(psm)$variables["peptide"] <- NA_character_
    expect_error(makeAdjacencyMatrix(psm))
    expect_error(adjacencyMatrix(psm))
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


test_that("ajacendyMatrix() accessor works", {
    psm <-  msdata::ident(full.names = TRUE, pattern = "TMT") |>
    PSM() |>
    filterPsmDecoy() |>
    filterPsmRank()
    ## set to binary to get the same result as the adjacencyMatrix,PSM
    ## accessor that created adj1
    cc <- ConnectedComponents(psm, binary = TRUE)
    ## not identical, a multiple PSMs per peptide
    adj1 <- adjacencyMatrix(psm)
    adj2 <- adjacencyMatrix(cc)
    expect_true(is(adj1, "sparseMatrix"))
    expect_true(is(adj2, "sparseMatrix"))
    expect_identical(sum(adj1@x), sum(adj2@x))
    expect_identical(colnames(adj1), colnames(adj2))
    ## remove duplicated sequences to make adjacency matrices
    ## identical
    psm <- psm[!duplicated(psm$sequence), ]
    cc <- ConnectedComponents(psm)
    adj1 <- adjacencyMatrix(psm)
    adj2 <- adjacencyMatrix(cc)
    expect_identical(adj1, adj2)
})


test_that("plotAjacendyMatrix() works", {
    cc <-  msdata::ident(full.names = TRUE, pattern = "TMT") |>
    PSM() |>
    filterPsmDecoy() |>
    filterPsmRank() |>
    ConnectedComponents()
    g <- plotAdjacencyMatrix(connectedComponents(cc, 672))
    expect_true(is(g, "igraph"))
    expect_identical(V(g)$color, c(NA_character_, "steelblue"))
    dev.off()
    g2 <- plotAdjacencyMatrix(connectedComponents(cc, 672), 1)
    expect_true(is(g2, "igraph"))
    expect_identical(V(g2)$color, c(0, 2))
    dev.off()
})

test_that("plotAjacendyMatrix() attributes work", {
    adj <-  msdata::ident(full.names = TRUE, pattern = "TMT") |>
    PSM() |>
    filterPsmDecoy() |>
    filterPsmRank() |>
    ConnectedComponents() |>
    connectedComponents(1082)
    ## ------------------------
    ## Test default colours
    g <- plotAdjacencyMatrix(adj)
    pep_nodes <- which(names(V(g)) %in% rownames(adj))
    prot_nodes <- which(names(V(g)) %in% colnames(adj))
    ## peptides are white/NA
    expect_true(all(is.na(V(g)$color[pep_nodes])))
    ## proteins are steelblue
    expect_true(all(V(g)$color[prot_nodes] == "steelblue"))
    ## ------------------------
    ## Test shapes
    expect_true(all(V(g)$shape[prot_nodes] == "square"))
    expect_true(all(V(g)$shape[pep_nodes] == "circle"))
    ## ---------------------------
    ## Test custom peptide colors
    exp_pep_cols <- c("red", "blue", "green", "red", "blue", "orange", "green")
    pep_cols <- c(rep("black", 10), exp_pep_cols)
    names(pep_cols) <- c(LETTERS[1:10], rownames(adj))
    pep_cols <- sample(pep_cols)
    g <- plotAdjacencyMatrix(adj, pepColors = pep_cols)
    ## proteins are still steelblue
    expect_true(all(V(g)$color[prot_nodes] == "steelblue"))
    ## peptides have set colours
    expect_identical(V(g)$color[pep_nodes], exp_pep_cols)
    ## ---------------------------
    ## Test custom protein colors
    exp_prot_cols <- c("red", "blue", "green", "orange")
    prot_cols <- c(rep("black", 10), exp_prot_cols)
    names(prot_cols) <- c(LETTERS[1:10], colnames(adj))
    prot_cols <- sample(prot_cols)
    g <- plotAdjacencyMatrix(adj, protColors = prot_cols)
    ## peptides are still white/NA
    expect_true(all(is.na(V(g)$color[pep_nodes])))
    ## proteins have set colours
    expect_identical(V(g)$color[prot_nodes], exp_prot_cols)
    ## ---------------------------------------
    ## Test custom protein and peptide colors
    g <- plotAdjacencyMatrix(adj, protColors = prot_cols,
                             pepColors = pep_cols)
    ## proteins have set colours
    expect_identical(V(g)$color[prot_nodes], exp_prot_cols)
    ## peptides have set colours
    expect_identical(V(g)$color[pep_nodes], exp_pep_cols)
})


test_that("ajacendyMatrix() with scores works", {
    psmdf <- data.frame(spectrum = c("sp1", "sp2", "sp3", "sp4", "sp5"),
                    sequence = c(c("A", LETTERS[1:4])),
                    protein = c("ProtA", "ProtA", "ProtB", "ProtC", "ProtD"),
                    decoy = rep(FALSE, 5),
                    rank = rep(1, 5),
                    score = 1:5)
    psm1 <- PSM(psmdf[-1, ], spectrum = "spectrum", peptide = "sequence",
                protein = "protein", decoy = "decoy", rank = "rank")
    psm2 <- PSM(psmdf, spectrum = "spectrum", peptide = "sequence",
                protein = "protein", decoy = "decoy", rank = "rank")
    psm3 <- PSM(psmdf[-1, ], spectrum = "spectrum", peptide = "sequence",
                protein = "protein", decoy = "decoy", rank = "rank", score = "score")
    psm4 <- PSM(psmdf, spectrum = "spectrum", peptide = "sequence",
                protein = "protein", decoy = "decoy", rank = "rank", score = "score")
    ## binary matrix
    adj1 <- makeAdjacencyMatrix(psm1)
    expect_identical(dimnames(adj1), list(c("A", "B", "C", "D"),
                                          c("ProtA", "ProtB", "ProtC", "ProtD")))
    expect_identical(sum(adj1), 4)
    expect_identical(adj1@x, rep(1, 4))
    ## A/Prot2 tallies 2 PSMs
    adj2 <- makeAdjacencyMatrix(psm2)
    expect_identical(dimnames(adj2), list(c("A", "B", "C", "D"),
                                          c("ProtA", "ProtB", "ProtC", "ProtD")))
    expect_identical(sum(adj2), 5)
    expect_identical(adj2@x, c(2, rep(1, 3)))
    expect_identical(adj1, makeAdjacencyMatrix(psm2, binary = TRUE))
    ## scores with single PSMs
    adj3 <- makeAdjacencyMatrix(psm3)
    expect_identical(dimnames(adj3), list(c("A", "B", "C", "D"),
                                          c("ProtA", "ProtB", "ProtC", "ProtD")))
    expect_identical(sum(adj3), 14)
    expect_identical(adj3@x, seq(2, 5, 1))
    expect_identical(adj1, makeAdjacencyMatrix(psm3, binary = TRUE))
    ## scores with one double PSMs
    adj4 <- makeAdjacencyMatrix(psm4)
    expect_identical(dimnames(adj4), list(c("A", "B", "C", "D"),
                                          c("ProtA", "ProtB", "ProtC", "ProtD")))
    expect_identical(sum(adj4), 15)
    expect_identical(adj4@x, c(3, seq(3, 5, 1)))
    expect_identical(adj1, makeAdjacencyMatrix(psm4, binary = TRUE))
})

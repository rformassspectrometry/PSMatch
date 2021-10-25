psmdf <- data.frame(spectrum = paste0("sp", 1:10),
                    sequence = replicate(10,
                                         paste(sample(getAminoAcids()[-1, "AA"], 10),
                                               collapse = "")),
                    protein = sample(paste0("Prot", LETTERS[1:7]), 10,
                                     replace = TRUE),
                    decoy = rep(FALSE, 10),
                    rank = rep(1, 10),
                    score = runif(10))

test_that("Test PSM construction from data.frame", {
    psm <- PSM(psmdf)
    expect_true(validObject(psm))
    expect_identical(nrow(psm), nrow(psmdf))
    expect_identical(ncol(psm), ncol(psmdf))
    expect_true(sum(is.na(psmVariables(psm))) == 5)
    psm2 <- PSM(psm, decoy = "decoy", rank = "rank",
                spectrum = "spectrum", peptide = "sequence",
                protein = "protein")
    expect_true(sum(is.na(psmVariables(psm2))) == 0)
})


test_that("Test PSM filtering", {
    psm <- PSM(psmdf)
    psm2 <- PSM(psm, decoy = "decoy", rank = "rank",
                spectrum = "spectrum", peptide = "sequence",
                protein = "protein")
    ## none of these filters on decoy should change the data
    expect_identical(psm, filterPsmDecoy(psm, decoy = NULL))
    expect_identical(psm, filterPsmDecoy(psm, decoy = "decoy"))
    expect_identical(psm2, filterPsmDecoy(psm2))
    expect_identical(psm2, filterPsmDecoy(psm2, decoy = "decoy"))
    expect_identical(psm2, filterPsmDecoy(psm2, decoy = NULL))
    ## none of these filters on rank should change the data
    expect_identical(psm, filterPsmRank(psm, rank = NULL))
    expect_identical(psm, filterPsmRank(psm, rank = "rank"))
    expect_identical(psm2, filterPsmRank(psm2))
    expect_identical(psm2, filterPsmRank(psm2, rank = NULL))
    expect_identical(psm2, filterPsmRank(psm2, rank = "rank"))
    ## none of these filters on uniqueness of hits should change the data
    expect_identical(psm, filterPsmNonProteotypic(psm, peptide = NULL))
    expect_identical(psm, filterPsmNonProteotypic(psm, protein = NULL))
    expect_identical(psm, filterPsmNonProteotypic(psm, peptide = "sequence",
                                                  protein = "protein"))
    expect_identical(psm2, filterPsmNonProteotypic(psm2))
    expect_identical(psm2, filterPsmNonProteotypic(psm2, peptide = NULL))
    expect_identical(psm2, filterPsmNonProteotypic(psm2, protein = NULL))
    expect_identical(psm2, filterPsmNonProteotypic(psm2, peptide = "sequence",
                                                   protein = "protein"))
})


test_that("Test PSM construction from msid", {
})

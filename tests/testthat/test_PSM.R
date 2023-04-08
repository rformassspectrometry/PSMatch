## Data for construction from mzid. Strangely, the test fails when
## creating the object with the mzID parser due to parallel
## processing.
f <- msdata::ident(full.names = TRUE, pattern = "TMT")
psm_mzR <- PSM(f, parser = "mzR")
psm_mzID <- PSM(f, parser = "mzID")

## No changes expected when filtering, given that there are no decoy
## peptides, all PSMs have rank 1, and all peptides match unique
## proteins.
psmdf0 <- data.frame(spectrum = paste0("sp", 1:10),
                     sequence = replicate(10,
                                          paste(sample(getAminoAcids()[-1, "AA"], 10),
                                                collapse = "")),
                     protein = paste0("Prot", LETTERS[1:10]),
                     decoy = rep(FALSE, 10),
                     rank = rep(1, 10),
                     score = runif(10),
                     fdr = runif(10))

## Filtering will have an effect here. There will be 12 PSMs left
## after filtering for decoys , 10 PSMs after filtering for rank, 5
## when keeping shared peptides. When filtering decoy and high rank
## PSMs, we also get rid of all the non-unique peptides, which leaves
## 10 PSMs.
set.seed(123) ## for score and fdr
psmdf1 <- data.frame(spectrum = paste0("sp", 1:15),
                     sequence = c(psmdf0$sequence,
                                  psmdf0$sequence[1:5]),
                     protein = paste0("Prot", LETTERS[1:15]),
                     decoy = c(rep(FALSE, 12), rep(TRUE, 3)),
                     rank = c(rep(1, 10), rep(2, 5)),
                     score = runif(15),
                     fdr = runif(15))
not_decoy <- sum(!psmdf1$decoy)   ## 12
rank_one <- sum(psmdf1$rank == 1) ## 10
unique_pep <- sum(table(psmdf1$sequence) == 1) ## 5
fdr_pep_1 <- sum(psmdf1$fdr < 0.05) ## 1
fdr_pep_5 <- sum(psmdf1$fdr < 0.5) ## 1

psmdf2 <- data.frame(spectrum = rep(paste0("sp", 1:5),
                                   c(2, 2, 1, 1, 1)),
                    sequence = psmdf0$sequence[1:7],
                    protein = c("ProtA", "ProtB",
                                "ProtA", "ProtA",
                                paste0("Prot", LETTERS[3:5])),
                    decoy = rep(FALSE, 7),
                    rank = rep(1, 7))

test_that("Test PSM construction from data.frame", {
    psm <- PSM(psmdf0)
    expect_true(validObject(psm))
    expect_identical(nrow(psm), nrow(psmdf0))
    expect_identical(ncol(psm), ncol(psmdf0))
    expect_true(sum(is.na(psmVariables(psm))) == 7)
    psm2 <- PSM(psm, decoy = "decoy", rank = "rank",
                spectrum = "spectrum", peptide = "sequence",
                protein = "protein")
    expect_true(sum(is.na(psmVariables(psm2))) == 2)
    psm2 <- PSM(psm, decoy = "decoy", rank = "rank",
                spectrum = "spectrum", peptide = "sequence",
                protein = "protein", score = "score",
                fdr = "fdr")
    expect_true(sum(is.na(psmVariables(psm2))) == 0)
    expect_error(psmVariables(psm, "error"))
    expect_error(psmVariables(psm2, "error"))
})

test_that("Test PSM filtering (no change)", {
    psm <- PSM(psmdf0)
    expect_true(validObject(psm))
    psm2 <- PSM(psm, decoy = "decoy", rank = "rank",
                spectrum = "spectrum", peptide = "sequence",
                protein = "protein")
    expect_true(validObject(psm2))
    ## none of these filters on decoy should change the data
    expect_identical(psm, filterPsmDecoy(psm, decoy = NA))
    expect_identical(psm, filterPsmDecoy(psm, decoy = NULL))
    expect_identical(psm, filterPsmDecoy(psm, decoy = "decoy"))
    expect_identical(psm2, filterPsmDecoy(psm2))
    expect_identical(psm2, filterPsmDecoy(psm2, decoy = "decoy"))
    expect_identical(psm2, filterPsmDecoy(psm2, decoy = NULL))
    ## none of these filters on rank should change the data
    expect_identical(psm, filterPsmRank(psm, rank = NA))
    expect_identical(psm, filterPsmRank(psm, rank = NULL))
    expect_identical(psm, filterPsmRank(psm, rank = "rank"))
    expect_identical(psm2, filterPsmRank(psm2))
    expect_identical(psm2, filterPsmRank(psm2, rank = NA))
    expect_identical(psm2, filterPsmRank(psm2, rank = NULL))
    expect_identical(psm2, filterPsmRank(psm2, rank = "rank"))
    ## none of these filters on uniqueness of hits should change the data
    expect_identical(psm, filterPsmShared(psm, peptide = NA))
    expect_identical(psm, filterPsmShared(psm, peptide = NULL))
    expect_identical(psm, filterPsmShared(psm, protein = NA))
    expect_identical(psm, filterPsmShared(psm, protein = NULL))
    expect_identical(psm, filterPsmShared(psm, peptide = "sequence",
                                          protein = "protein"))
    expect_identical(psm2, filterPsmShared(psm2))
    expect_identical(psm2, filterPsmShared(psm2, peptide = NA))
    expect_identical(psm2, filterPsmShared(psm2, peptide = NULL))
    expect_identical(psm2, filterPsmShared(psm2, protein = NA))
    expect_identical(psm2, filterPsmShared(psm2, protein = NULL))
    expect_identical(psm2, filterPsmShared(psm2, peptide = "sequence",
                                           protein = "protein"))
})

test_that("Test PSM construction from mzid files", {
    ## Valid objects
    expect_true(validObject(psm_mzR))
    expect_true(validObject(psm_mzID))
    ## No missing PSM variables
    expect_equal(sum(is.na(psmVariables(psm_mzR))), 2)
    expect_equal(sum(is.na(psmVariables(psm_mzID))), 2)
    nms <- c("spectrum", "peptide", "protein", "decoy", "rank")
    expect_equal(sum(is.na(psmVariables(psm_mzR)[nms])), 0)
    expect_equal(sum(is.na(psmVariables(psm_mzID)[nms])), 0)
    expect_true(is.na(psmVariables(psm_mzR)["score"]))
    expect_true(is.na(psmVariables(psm_mzID)["score"]))
    expect_true(is.na(psmVariables(psm_mzR)["fdr"]))
    expect_true(is.na(psmVariables(psm_mzID)["fdr"]))
    expect_identical(names(psmVariables(psm_mzID)),
                     names(psmVariables(psm_mzR)))
})


test_that("Test PSM show", {
    expect_null(show(psm_mzR))
    expect_null(show(psm_mzID))
    expect_null(show(PSM(psmdf0)))
    expect_null(show(PSM(psmdf1)))
    ## show object with ncol() <= 4
    expect_null(show(PSM(psmdf0[, 1:3])))
})

test_that("psmVariables accessor works", {
    expect_identical(psmVariables(psm_mzR)["peptide"],
                     psmVariables(psm_mzR, "peptide"))
})

test_that("Test PSM filtering", {
    psm <- PSM(psmdf1, decoy = "decoy", rank = "rank",
               spectrum = "spectrum", peptide = "sequence",
               protein = "protein", fdr = "fdr")
    expect_true(validObject(psm))
    expect_identical(nrow(filterPsmDecoy(psm)), not_decoy)
    expect_identical(nrow(filterPsmRank(psm)), rank_one)
    expect_identical(nrow(filterPsmShared(psm)), unique_pep)
    expect_identical(nrow(filterPSMs(psm)), 10L)
    expect_identical(nrow(filterPsmFdr(psm)), fdr_pep_1)
    expect_identical(nrow(filterPsmFdr(psm, FDR = 0.5)), fdr_pep_5)
})


test_that("reducePSMs works", {
    psm <- PSM(psmdf2)
    expect_identical(nrow(psm), 7L)
    expect_true(is.na(reduced(psm)))
    reduced(psm) <- FALSE
    expect_false(reduced(psm))
    rpsm <- reducePSMs(psm, psm[["spectrum"]])
    expect_null(show(rpsm))
    expect_true(reduced(rpsm))
    expect_identical(nrow(rpsm), 5L)
    expect_identical(rpsm[["spectrum"]], unique(psm[["spectrum"]]))
    expect_identical(rpsm[["decoy"]], rep(FALSE, 5))
    expect_identical(rpsm[["rank"]], rep(1, 5))
    expect_equal(rpsm[["protein"]],
                 IRanges::CharacterList(list(sp1 = psmdf2$protein[1:2],
                                             sp2 = psmdf2$protein[3],
                                             sp3 = psmdf2$protein[5],
                                             sp4 = psmdf2$protein[6],
                                             sp5 = psmdf2$protein[7])))
    expect_equal(rpsm[["sequence"]],
                 IRanges::CharacterList(list(sp1 = psmdf2$sequence[1:2],
                                             sp2 = psmdf2$sequence[3:4],
                                             sp3 = psmdf2$sequence[5],
                                             sp4 = psmdf2$sequence[6],
                                             sp5 = psmdf2$sequence[7])))
})

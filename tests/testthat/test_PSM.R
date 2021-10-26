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
                     score = runif(10))

## Filtering will have an effect here. There will be 12 PSMs left
## after filtering for decoys , 10 PSMs after filtering for rank, 5
## when keeping proteotypic peptides. When filtering decoy and high
## rank PSMs, we also get rid of all the non-unique peptides, which
## leaves 10 PSMs.
psmdf1 <- data.frame(spectrum = paste0("sp", 1:15),
                     sequence = c(psmdf0$sequence,
                                  psmdf0$sequence[1:5]),
                     protein = paste0("Prot", LETTERS[1:15]),
                     decoy = c(rep(FALSE, 12), rep(TRUE, 3)),
                     rank = c(rep(1, 10), rep(2, 5)),
                     score = runif(15))
not_decoy <- sum(!psmdf1$decoy)   ## 12
rank_one <- sum(psmdf1$rank == 1) ## 10
unique_pep <- sum(table(psmdf1$sequence) == 1) ## 5


test_that("Test PSM construction from data.frame", {
    psm <- PSM(psmdf0)
    expect_true(validObject(psm))
    expect_identical(nrow(psm), nrow(psmdf0))
    expect_identical(ncol(psm), ncol(psmdf0))
    expect_true(sum(is.na(psmVariables(psm))) == 5)
    psm2 <- PSM(psm, decoy = "decoy", rank = "rank",
                spectrum = "spectrum", peptide = "sequence",
                protein = "protein")
    expect_true(sum(is.na(psmVariables(psm2))) == 0)
})

test_that("Test PSM filtering (no change)", {
    psm <- PSM(psmdf0)
    expect_true(validObject(psm))
    psm2 <- PSM(psm, decoy = "decoy", rank = "rank",
                spectrum = "spectrum", peptide = "sequence",
                protein = "protein")
    expect_true(validObject(psm2))
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

test_that("Test PSM filtering", {
    psm <- PSM(psmdf1, decoy = "decoy", rank = "rank",
               spectrum = "spectrum", peptide = "sequence",
               protein = "protein")
    expect_true(validObject(psm))
    expect_identical(nrow(filterPsmDecoy(psm)), not_decoy)
    expect_identical(nrow(filterPsmRank(psm)), rank_one)
    expect_identical(nrow(filterPsmNonProteotypic(psm)), unique_pep)
    expect_identical(nrow(filterPSMs(psm)), 10L)
})

test_that("Test PSM construction from mzid", {
    ## Valid objects
    expect_true(validObject(psm_mzR))
    expect_true(validObject(psm_mzID))
    ## No missing PSM variables
    expect_equal(sum(is.na(psmVariables(psm_mzR))), 0)
    expect_equal(sum(is.na(psmVariables(psm_mzID))), 0)
    expect_identical(names(psmVariables(psm_mzID)),
                     names(psmVariables(psm_mzR)))
})

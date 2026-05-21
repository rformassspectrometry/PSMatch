library("Spectra")

##  Helpers

.makeSpectra <- function(mz, intensity, sequence = NULL,
                         precursorMz = NA_real_) {
    df <- DataFrame(msLevel = 2L, rtime = 100,
                    precursorMz = precursorMz)
    df$mz <- list(mz)
    df$intensity <- list(intensity)
    if (!is.null(sequence))
        df$sequence <- sequence
    Spectra(df)
}

##  checkABpresence

test_that("checkABpresence() returns TRUE when a2 and b2 are both present", {
    seq <- "PEPTIDE"
    frags <- calculateFragments(seq, type = c("a", "b"))
    sp <- .makeSpectra(sort(frags$mz), rep(1e4, nrow(frags)), seq)
    expect_true(checkABpresence(sp))
})

test_that("checkABpresence() returns FALSE when a2 is absent", {
    seq <- "PEPTIDE"
    frags <- calculateFragments(seq, type = "b")
    sp <- .makeSpectra(sort(frags$mz), rep(1e4, nrow(frags)), seq)
    expect_false(checkABpresence(sp))
})

test_that("checkABpresence() works with pre-computed fragments", {
    sp <- .makeSpectra(c(100, 200, 300), rep(1e4, 3))
    labs_ab  <- list(c("a1", "a2", "b1", "b2", "b3"))
    labs_b   <- list(c("b1", "b2", "b3"))
    expect_true(checkABpresence(sp, fragments = labs_ab))
    expect_false(checkABpresence(sp, fragments = labs_b))
})

test_that("checkABpresence() handles multiple spectra", {
    df <- DataFrame(msLevel = c(2L, 2L), rtime = c(100, 101))
    df$mz <- list(c(100, 200), c(100, 200))
    df$intensity <- list(c(1e4, 1e4), c(1e4, 1e4))
    sp2 <- Spectra(df)
    labs <- list(c("a2", "b2"), c("b2", "b3"))
    expect_equal(checkABpresence(sp2, fragments = labs), c(TRUE, FALSE))
})

##  checkXYpresence

test_that("checkXYpresence() works with pre-computed fragments", {
    sp <- .makeSpectra(c(100, 200, 300), rep(1e4, 3))
    labs_xy <- list(c("x2", "y2", "y3"))
    labs_y  <- list(c("y2", "y3", "y4"))
    expect_true(checkXYpresence(sp, fragments = labs_xy))
    expect_false(checkXYpresence(sp, fragments = labs_y))
})

##  checkParentIonIntensity

test_that("checkParentIonIntensity() returns 1 when precursor = base peak", {
    sp <- .makeSpectra(c(100, 200, 300, 400, 500),
                       c(100, 100, 100, 100, 1e6),
                       precursorMz = 500.0)
    expect_equal(checkParentIonIntensity(sp, tolerance = 0), 1)
})

test_that("checkParentIonIntensity() returns 0 when precursor is absent", {
    sp <- .makeSpectra(c(100, 200, 300, 400, 500),
                       c(100, 100, 100, 100, 1e6),
                       precursorMz = 999.0)
    expect_equal(checkParentIonIntensity(sp, tolerance = 0), 0)
})

##  checkOverlap

test_that("checkOverlap() returns TRUE when b/y-ions fully overlap", {
    seq <- "PEPTIDE"
    frags <- calculateFragments(seq, type = c("b", "y"))
    sp <- .makeSpectra(sort(frags$mz), rep(1e4, nrow(frags)), seq)
    labs <- suppressWarnings(labelFragments(sp, type = c("b", "y")))
    expect_true(checkOverlap(sp, fragments = labs, strippedSeq = seq))
})

test_that("checkOverlap() returns FALSE when coverage gap is present", {
    seq <- "PEPTIDE"
    frags <- calculateFragments(seq, type = c("b", "y"))
    frags_sparse <- frags[frags$ion %in% c("b1", "b2", "y1", "y2"), ]
    sp <- .makeSpectra(sort(frags_sparse$mz),
                       rep(1e4, nrow(frags_sparse)), seq)
    labs <- suppressWarnings(labelFragments(sp, type = c("b", "y")))
    expect_false(checkOverlap(sp, fragments = labs, strippedSeq = seq))
})

##  checkShiftConsistency

test_that("checkShiftConsistency() returns NA for unmodified sequences", {
    seq <- "PEPTIDE"
    frags <- calculateFragments(seq)
    sp <- .makeSpectra(sort(frags$mz), rep(1e4, nrow(frags)), seq)
    labs <- suppressWarnings(labelFragments(sp))
    result <- checkShiftConsistency(sp, fragments = labs)
    expect_true(is.na(result))
})

##  validatePSM

test_that("validatePSM() returns a data.frame with expected columns", {
    data("psmBoekweg")
    data("spBoekweg")
    psmBoekweg$pkey <- paste0(
        basename(psmBoekweg$filename),
        sub("^.+scan=", "::", psmBoekweg$scannr))
    sp <- Spectra::joinSpectraData(spBoekweg, psmBoekweg, by.x = "pkey")
    seq_var <- psmVariables(psmBoekweg)[["peptide"]]
    fdr_var <- psmVariables(psmBoekweg)[["fdr"]]
    res <- suppressWarnings(
        validatePSM(sp[9:10], peptideVariable = seq_var,
                    fdr = fdr_var))
    expect_s3_class(res, "data.frame")
    expect_true(all(c("spectrumId", "scanIndex", "peptide", "canonicalSeq",
                       "fdr", "a2b2", "x2y2", "byOverlap",
                       "shiftConsistency", "parentIonInt",
                       "precursorPurity") %in% names(res)))
    ## Only identified PSMs should appear
    expect_true(all(!is.na(res$peptide)))
})
#' # each base::graphics plot function must be wrapped by an anonymous function
#' # that could be called by `vdiffr::expect_doppelganger()`
#' Run devtools::test_active_file(file = "tests/testthat/test_plotSpectraPTM.R")

library("Spectra")

sp <- DataFrame(
    msLevel = 2L,
    rtime = 2345,
    sequence = "HIGFEGDSIGR",
    dataOrigin = "testfile.mzML",
    scanIndex = 1L,
    charge = 2L
)
sp$mz <- list(c(
    223.1583, 251.15432, 308.168017, 455.24801, 604.30949,
    641.30842, 667.2244, 778.30164, 813.34935, 923.350391,
    995.45281, 1017.43394, 1065.46197, 1112.5069, 1130.5874
))
sp$intensity <- list(c(
    83000, 65000, 190000, 379000, 281000, 112000, 39000,
    139000, 1015000, 63000, 58000, 1960000, 240000,
    1338000, 40700
))
spectra <- Spectra(sp)

test_that("plotSpectraPTM works with deltaMz = TRUE", {
    expect_doppelganger(
        "deltaMz-true",
        function() {
            plotSpectraPTM(
                spectra,
                type = c("a", "b", "c", "x", "y", "z"),
                deltaMz = TRUE,
                z = 1
            )
        }
    )
})

test_that("plotSpectraPTM works with deltaMz = FALSE", {
    expect_doppelganger(
        "deltaMz-false",
        function() {
            plotSpectraPTM(
                spectra,
                type = c("a", "b", "c", "x", "y", "z"),
                deltaMz = FALSE,
                z = 1
            )
        }
    )
})

# test_that("plotSpectraPTM works with variable modifications", {
#     expect_doppelganger(
#         "one-ptm-deltaMz-true",
#         function() {
#             plotSpectraPTM(
#                 spectra,
#                 type = c("a", "b", "c", "x", "y", "z"),
#                 variable_modifications = c(S = 79.996),
#                 deltaMz = TRUE
#             )
#         }
#     )
# })

test_that("plotSpectraPTM works with different col", {
    # We're using fixed colors here for reproducibility
    expect_doppelganger(
        "diff-col",
        function() {
            plotSpectraPTM(
                spectra,
                col = c(y = "red", b = "blue", acxy = "orange", other = "violet"),
                type = c("a", "b", "c", "x", "y", "z"),
                deltaMz = FALSE,
                z = 1
            )
        }
    )
})

test_that("plotSpectraPTM expands per-spectrum z list in sync with variableModifications", {
    ## Two spectra with different precursor charges.
    ## variableModifications = c(A = 1.0, Q = 1.0) expands each by one combination:
    ##   "ACE" -> "ACE", "A[+1.0]CE"   (2 spectra)
    ##   "PQR" -> "PQR", "PQ[+1.0]R"   (2 spectra)
    ## z = list(1:2, 1:3) must be replicated to list(1:2, 1:2, 1:3, 1:3).
    ## Without the replication, labelFragments throws "subscript out of bounds".
    sp2 <- DataFrame(
        msLevel = c(2L, 2L),
        rtime = c(100, 200),
        sequence = c("ACE", "PQR"),
        dataOrigin = c("f.mzML", "f.mzML"),
        scanIndex = c(1L, 2L),
        charge = c(2L, 3L)
    )
    sp2$mz <- list(c(100.0, 200.0, 300.0), c(150.0, 250.0, 350.0))
    sp2$intensity <- list(c(1000.0, 2000.0, 3000.0), c(1500.0, 2500.0, 3500.0))
    sps2 <- Spectra(sp2)

    expect_no_error(
        plotSpectraPTM(sps2,
            z = list(1:2, 1:3),
            variableModifications = c(A = 1.0, Q = 1.0),
            addCarbamidomethyl = FALSE,
            deltaMz = FALSE
        )
    )
})

test_that("plotSpectraPTM works with USI = FALSE", {
    expect_doppelganger(
        "USI-true",
        function() {
            plotSpectraPTM(
                spectra,
                type = c("a", "b", "c", "x", "y", "z"),
                USI = FALSE,
                z = 1
            )
        }
    )
})

test_that("plotSpectraPTM works with allCharges = FALSE", {
    expect_doppelganger(
        "allCharges-false",
        function() {
            plotSpectraPTM(
                spectra,
                type = c("a", "b", "c", "x", "y", "z"),
                allCharges = FALSE,
                deltaMz = FALSE
            )
        }
    )
})

test_that("plotSpectraPTM works with custom z parameter", {
    sp2 <- c(spectra, spectra)
    sp2$precursorCharge <- c(2L, 3L)
    expect_doppelganger(
        "custom-z",
        function() {
            plotSpectraPTM(
                sp2,
                type = c("a", "b", "c", "x", "y", "z"),
                z = c(2, 3),
                deltaMz = FALSE
            )
        }
    )
})

test_that("plotSpectraPTM errors with z and variableModifications", {
    expect_error(
        plotSpectraPTM(
            spectra,
            z = c(2),
            variableModifications = c(H = 15.994915)
        ),
        "Cannot use both 'z' and 'variableModifications'"
    )
})

test_that("plotSpectraPTM errors with wrong length z", {
    sp2 <- c(spectra, spectra)
    expect_error(
        plotSpectraPTM(sp2, z = c(2)),
        "'z' must be NULL or a numeric vector of length equal to length\\(x\\)"
    )
})
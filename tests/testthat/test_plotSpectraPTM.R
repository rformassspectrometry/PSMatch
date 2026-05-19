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
sp$mz <- list(c(223.1583, 251.15432, 308.168017, 455.24801, 604.30949,
                641.30842, 667.2244, 778.30164, 813.34935, 923.350391,
                995.45281, 1017.43394, 1065.46197, 1112.5069, 1130.5874))
sp$intensity <- list(c(83000, 65000, 190000, 379000, 281000, 112000, 39000,
                       139000, 1015000, 63000, 58000, 1960000, 240000,
                       1338000, 40700))
spectra <- Spectra(sp)

test_that("plotSpectraPTM works with deltaMz = TRUE", {
    expect_doppelganger(
        "deltaMz-true",
        function() {
            plotSpectraPTM(
                spectra,
                type = c("a", "b", "c", "x", "y", "z"),
                deltaMz = TRUE
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
                deltaMz = FALSE
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
                deltaMz = FALSE
            )
        }
    )
})

test_that("plotSpectraPTM works with USI = FALSE", {
    expect_doppelganger(
        "USI-true",
        function() {
            plotSpectraPTM(
                spectra,
                type = c("a", "b", "c", "x", "y", "z"),
                USI = FALSE
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
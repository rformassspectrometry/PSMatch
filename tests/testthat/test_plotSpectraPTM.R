#' # each base::graphics plot function must be wrapped by an anonymous function
#' # that could be called by `vdiffr::expect_doppelganger()`
#' Run devtools::test_active_file(file = "tests/testthat/test_plotSpectraPTM.R")

test_that("plotSpectraPTM works with deltaMz = TRUE", {
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

    # We're using fixed colors here for reproducibility
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

    # We're using fixed colors here for reproducibility
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
#     sp <- DataFrame(
#         msLevel = 2L,
#         rtime = 2345,
#         sequence = "HIGFEGDSIGR",
#         dataOrigin = "testfile.mzML",
#         scanIndex = 1L,
#         charge = 2L
#     )
#     sp$mz <- list(c(223.1583, 251.15432, 308.168017, 455.24801, 604.30949,
#                     641.30842, 667.2244, 778.30164, 813.34935, 923.350391,
#                     995.45281, 1017.43394, 1065.46197, 1112.5069, 1130.5874))
#     sp$intensity <- list(c(83000, 65000, 190000, 379000, 281000, 112000, 39000,
#                            139000, 1015000, 63000, 58000, 1960000, 240000,
#                            1338000, 40700))
#     spectra <- Spectra(sp)
# 
#     # We're using fixed colors here for reproducibility
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
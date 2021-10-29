test_that("addFragments() works", {
    library("Spectra")
    seq <- "PQR"
    frags <- calculateFragments(seq)
    o <- order(frags$mz)
    sp <- DataFrame(msLevel = 2L, rtime = 2345, sequence = seq)
    sp$mz <- list(frags$mz[o])
    sp$intensity <- list(rep(1, 7))
    sp <- Spectra(sp)
    ## all fragments
    ans <- addFragments(sp)
    exp <- frags$ion[o]
    expect_identical(ans, exp)
    ## 2nd fragment missing
    sp$mz[[1]][2] <- sp$mz[[1]][2] * 1.1
    exp[2] <- NA
    ans <- addFragments(sp)
    expect_identical(ans, exp)
    ## 7th fragment missing
    sp$mz[[1]][7] <- sp$mz[[1]][7] * 1.1
    exp[7] <- NA
    ans <- addFragments(sp)
    expect_identical(ans, exp)
})

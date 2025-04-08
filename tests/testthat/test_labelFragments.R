test_that("labelFragments() works", {
    library("Spectra")
    seq <- "PQR"
    frags <- calculateFragments(seq)
    o <- order(frags$mz)
    sp <- DataFrame(msLevel = 2L, rtime = 2345, sequence = seq)
    sp$mz <- list(frags$mz[o])
    sp$intensity <- list(rep(1, 7))
    sp <- Spectra(sp)
    ## all fragments
    ans <- labelFragments(sp)[[1]]
    exp <- frags$ion[o]
    ## Remove attribute
    expect_identical(ans[1:length(ans)], exp)
    ## 2nd fragment missing
    sp$mz[[1]][2] <- sp$mz[[1]][2] * 1.1
    exp[2] <- NA
    ans <- labelFragments(sp)[[1]]
    expect_identical(ans[1:length(ans)], exp)
    ## 7th fragment missing
    sp$mz[[1]][7] <- sp$mz[[1]][7] * 1.1
    exp[7] <- NA
    ans <- labelFragments(sp)[[1]]
    expect_identical(ans[1:length(ans)], exp)
})


test_that("labelFragments() works with multiple Spectra", {
    library("Spectra")
    seq <- c("PQR", "ACE")
    frags_pqr <- calculateFragments(seq)[1:7,]
    frags_ace <- calculateFragments(seq)[8:13,]
    o_pqr <- order(frags_pqr$mz)
    o_ace <- order(frags_ace$mz)
    sp <- DataFrame(msLevel = c(2L, 2L), rtime = c(2345, 2346), sequence = seq)
    sp$mz <- c(list(frags_pqr$mz[o_pqr]), list(frags_ace$mz[o_ace]))
    sp$intensity <- c(list(rep(1, 7)), list(rep(1, 6)))
    sp <- Spectra(sp)
    ## all fragments
    ans <- labelFragments(sp)
    ## Number of elements equal the possibilities of peptide sequences
    ## This instance: no mod, 1 mod on Q: 2 possibilities
    expect_equal(length(ans), 2)
    expect_identical(names(ans), seq)
})

test_that("labelFragments() works with modifications", {
    library("Spectra")
    seq <- "PQR"
    frags <- calculateFragments(seq)
    o <- order(frags$mz)
    sp <- DataFrame(msLevel = 2L, rtime = 2345, sequence = seq)
    sp$mz <- list(frags$mz[o])
    sp$intensity <- list(rep(1, 7))
    sp <- Spectra(sp)
    ## all fragments
    ans <- labelFragments(sp, variable_modifications = c(Q = 45))
    ## Number of elements equal the possibilities of peptide sequences
    ## This instance: no mod, 1 mod on Q: 2 possibilities
    expect_equal(length(ans), 2)
})


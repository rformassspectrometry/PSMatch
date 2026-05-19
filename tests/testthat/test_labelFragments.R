library("Spectra")

test_that("labelFragments() works", {
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
    seq <- c("PQR", "ACE")
    frags_pqr <- calculateFragments(seq, addCarbamidomethyl = FALSE)[1:7,]
    frags_ace <- calculateFragments(seq, addCarbamidomethyl = FALSE)[8:13,]
    o_pqr <- order(frags_pqr$mz)
    o_ace <- order(frags_ace$mz)
    sp <- DataFrame(msLevel = c(2L, 2L), rtime = c(2345, 2346), sequence = seq)
    sp$mz <- c(list(frags_pqr$mz[o_pqr]), list(frags_ace$mz[o_ace]))
    sp$intensity <- c(list(rep(1, 7)), list(rep(1, 6)))
    sp <- Spectra(sp)
    ## all fragments
    ans <- labelFragments(sp, addCarbamidomethyl = FALSE)
    ## Number of elements equal the possibilities of peptide sequences
    ## This instance: no mod, 1 mod on Q: 2 possibilities
    expect_equal(length(ans), 2)
    expect_identical(names(ans), seq)
})

test_that("labelFragments() works with modifications", {
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

test_that("labelFragments() works with what = 'mz'", {
    seq <- "PQR"
    frags <- calculateFragments(seq)
    o <- order(frags$mz)
    sp <- DataFrame(msLevel = 2L, rtime = 2345, sequence = seq)
    sp$mz <- list(frags$mz[o])
    sp$intensity <- list(rep(1, 7))
    sp <- Spectra(sp)
    ## all fragments
    ans_mz<- labelFragments(sp, what = "mz")
    expect_identical(frags[o,"mz"], ans_mz[[1]])
    ans_ions <- labelFragments(sp, what = "ion")
    expect_identical(frags[o,"ion"], ans_ions[[1]])
})

test_that("labelFragments() works with custom z values", {
    seq <- c("PQR", "ACE")
    ## Create fragments with different charge states
    frags_pqr <- calculateFragments(seq[1], z = 1:2,
                                     addCarbamidomethyl = FALSE)
    frags_ace <- calculateFragments(seq[2], z = 1,
                                     addCarbamidomethyl = FALSE)
    ## Order fragments by mz
    o_pqr <- order(frags_pqr$mz)
    o_ace <- order(frags_ace$mz)
    ## Create spectra
    sp <- DataFrame(msLevel = c(2L, 2L), rtime = c(2345, 2346),
                    sequence = seq, precursorCharge = c(2L, 1L))
    sp$mz <- list(frags_pqr$mz[o_pqr][1:10],
                  frags_ace$mz[o_ace][1:5])
    sp$intensity <- list(rep(1, 10), rep(1, 5))
    sp <- Spectra(sp)
    ## Test with custom z values (numeric vector)
    ans <- labelFragments(sp, z = c(2, 1), allCharges = TRUE,
                          addCarbamidomethyl = FALSE)
    expect_equal(length(ans), 2)
    expect_equal(names(ans), seq)
})

test_that("labelFragments() works with allCharges = FALSE", {
    seq <- "PQR"
    ## Create fragments with charge 1 only
    frags <- calculateFragments(seq, z = 1,
                                addCarbamidomethyl = FALSE)
    o <- order(frags$mz)
    sp <- DataFrame(msLevel = 2L, rtime = 2345, sequence = seq,
                    precursorCharge = 2L)
    sp$mz <- list(frags$mz[o])
    sp$intensity <- list(rep(1, length(o)))
    sp <- Spectra(sp)
    ## Test with allCharges = FALSE (only charge 1)
    ans <- labelFragments(sp, z = c(1), allCharges = FALSE,
                          addCarbamidomethyl = FALSE)
    expect_equal(length(ans), 1)
})
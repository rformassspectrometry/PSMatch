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

test_that("labelFragments() respects per-spectrum z list", {
    ## Build reference fragments without neutral losses for deterministic mz values
    frags_ace <- calculateFragments("ACE", z = 1:2, addCarbamidomethyl = FALSE,
                                    neutralLoss = NULL)
    frags_pqr <- calculateFragments("PQR", z = 1:3, addCarbamidomethyl = FALSE,
                                    neutralLoss = NULL)

    ## Pick one ion unique to each higher charge state
    z2_row <- frags_ace[frags_ace$z == 2, ][1L, ]
    z3_row <- frags_pqr[frags_pqr$z == 3, ][1L, ]

    sp <- DataFrame(
        msLevel   = c(2L, 2L),
        rtime     = c(100, 200),
        sequence  = c("ACE", "PQR")
    )
    sp$mz       <- list(z2_row$mz, z3_row$mz)
    sp$intensity <- list(1.0, 1.0)
    sp <- Spectra(sp)

    ## Per-spectrum z list: both higher-charge ions are matched
    ans <- labelFragments(sp, z = list(1:2, 1:3), addCarbamidomethyl = FALSE,
                          neutralLoss = NULL)
    expect_identical(ans[[1L]], z2_row$ion)
    expect_identical(ans[[2L]], z3_row$ion)

    ## Scalar z = 1: those mz values fall outside the z=1 search space
    ans_z1 <- labelFragments(sp, z = 1L, addCarbamidomethyl = FALSE,
                             neutralLoss = NULL)
    expect_true(all(is.na(ans_z1[[1L]])))
    expect_true(all(is.na(ans_z1[[2L]])))
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
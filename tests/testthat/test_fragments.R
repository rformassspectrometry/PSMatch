test_that("calculateFragments", {
    pqr <- data.frame(
        mz = c(70.065, 198.124,   # a
               98.060, 226.119,   # b
               115.087, 243.145,  # c
               201.098, 329.157,  # x
               175.119, 303.178,  # y
               158.092, 286.151,  # z
               157.108, 285.167,  # y_
               286.151),          # y*
        ion = c(paste0(rep(c("a", "b", "c", "x", "y", "z"), each=2),
                       rep(1:2, times = 6)),
                paste0("y", 1:2, "_"), "y2*"),
        type = c(rep(c("a", "b", "c", "x", "y", "z", "y_"), each=2),
                 "y*"),
        pos = c(rep(1:2, 7), 2),
        z = 1,
        seq = c(rep(c("P", "PQ"), 3),
                rep(c("R", "QR"), 4),
                "QR"),
        peptide = rep("PQR", 15),
        stringsAsFactors=FALSE)
    
    ace <- data.frame(
        mz = c(22.528, 102.544,  # a
               36.526, 116.541,  # b
               45.039, 125.054,  # c
               87.523, 167.539,  # x
               74.534, 154.549,  # y
               66.021, 146.036), # z
        ion = paste0(rep(c("a", "b", "c", "x", "y", "z"), each=2),
                     rep(1:2, times = 6)),
        type = rep(c("a", "b", "c", "x", "y", "z"), each=2),
        pos = rep(1:2, 6),
        z = 2,
        seq = c(rep(c("A", "AC"), 3),
                rep(c("E", "CE"), 3)),
        peptide = rep("ACE", 12),
        stringsAsFactors=FALSE)
    
    expect_warning(calculateFragments("PQR", modifications = c(P=2)),
                   "'modifications' is deprecated, please use 'fixed_modifications' instead.")
    
    expect_equal(pqr[1:12,],
                 calculateFragments("PQR",
                                    type = c("a", "b", "c", "x", "y", "z"),
                                    neutralLoss = NULL, verbose = FALSE),
                 tolerance=1e-5)
    expect_equal(pqr[1:4,],
                 calculateFragments("PQR", type = c("a", "b"),
                                    neutralLoss = NULL, verbose = FALSE),
                 tolerance=1e-5)
    ## rownames always differ
    expect_equal(pqr[c(7:8, 11:12),],
                 calculateFragments("PQR", type = c("x", "z"),
                                    neutralLoss = NULL, verbose = FALSE),
                 check.attributes = FALSE, tolerance = 1e-5)
    
    ## neutral loss
    ## rownames always differ
    expect_equal(pqr[c(3:4, 9:10, 13:15),],
                 calculateFragments("PQR", verbose = FALSE),
                 check.attributes = FALSE, tolerance = 1e-5)
    
    ## neutral loss (water=cterm disabled),
    ## rownames always differ
    expect_equal(pqr[c(3:4, 9:10, 15),],
                 calculateFragments("PQR",
                                    neutralLoss = defaultNeutralLoss(disableWaterLoss = "Cterm"),
                                    verbose = FALSE),
                 check.attributes = FALSE, tolerance = 1e-5)
    
    ## neutral loss (ammonia=Q disabled),
    ## rownames always differ
    expect_equal(pqr[c(3:4, 9:10, 13:14),],
                 calculateFragments("PQR",
                                    neutralLoss = defaultNeutralLoss(disableAmmoniaLoss = "Q"),
                                    verbose = FALSE),
                 check.attributes = FALSE, tolerance = 1e-5)
    
    ## neutral loss + nterm mod, rownames always differ
    tpqr <- pqr[c(3:4, 9:10, 13:15),]
    tpqr$mz[1:2] <- tpqr$mz[1:2] + 229
    expect_equal(tpqr,
                 calculateFragments("PQR", fixed_modifications = c(C = 57.02146, Nterm = 229),
                                    verbose = FALSE),
                 check.attributes = FALSE, tolerance = 1e-5)
    
    ## neutral loss + nterm + cterm mod, rownames always differ
    tpqr$mz[3:7] <- tpqr$mz[3:7] - 100
    expect_equal(tpqr,
                 calculateFragments("PQR", fixed_modifications = c(C = 57.02146,
                                                             Nterm = 229,
                                                             Cterm = -100),
                                    verbose = FALSE),
                 check.attributes = FALSE, tolerance = 1e-5)
    
    expect_equal(ace,
                 calculateFragments("ACE", type = c("a", "b", "c", "x", "y", "z"),
                                    z = 2, neutralLoss = NULL, verbose = FALSE),
                 tolerance = 1e-5)
    expect_equal(ace[1:6,],
                 calculateFragments("ACE", type = letters[1:3], z = 2, verbose = FALSE),
                 tolerance = 1e-5)
    
    expect_error(calculateFragments("A"), "two or more residues")
    
    ## issue #200 (mz are not calculated correctly for terminal fixed_modifications
    ## and z > 1)
    p <- getAtomicMass()["p"]
    expect_equal(calculateFragments("AA", z = 2,
                                    fixed_modifications = c(Nterm = 10),
                                    type = "b")$mz - p,
                 (calculateFragments("AA", z = 1,
                                     fixed_modifications = c(Nterm = 10),
                                     type = "b")$mz - p )/ 2)
    expect_equal(calculateFragments("AA", z = 2, neutralLoss = NULL,
                                    fixed_modifications = c(Cterm = 10),
                                    type = "y")$mz - p,
                 (calculateFragments("AA", z = 1, neutralLoss = NULL,
                                     fixed_modifications = c(Cterm = 10),
                                     type = "y")$mz - p) / 2)
    
    ## See issue 573 in MSnbase (charge is ignored in neutral loss
    ## calculation)
    expect_equal(
        subset(calculateFragments("PEPTIDEE", z = 3, type = "b"), pos == 7L),
        data.frame(
            mz = c(261.4570693, 255.4535476),
            ion = c("b7", "b7_"),
            type = c("b", "b_"),
            pos = 7L,
            z = 3,
            seq = "PEPTIDE",
            peptide = "PEPTIDEE"
        ),
        check.attributes = FALSE # row.names differ
    )
    
})

test_that("defaultNeutralLoss", {
    expect_equal(defaultNeutralLoss(),
                 list(water = c("Cterm", "D", "E", "S", "T"),
                      ammonia = c("K", "N", "Q", "R")))
    expect_equal(defaultNeutralLoss(disableWaterLoss = c("T", "E", "S", "D")),
                 list(water = c("Cterm"), ammonia = c("K", "N", "Q", "R")))
    expect_equal(defaultNeutralLoss(disableWaterLoss = c("T", "E", "S", "D"),
                                    disableAmmoniaLoss = c("K", "Q")),
                 list(water = c("Cterm"), ammonia = c("N", "R")))
    expect_equal(defaultNeutralLoss(disableWaterLoss = c("Cterm",
                                                         "T", "E", "S", "D"),
                                    disableAmmoniaLoss = c("K", "N", "Q", "R")),
                 list(water = character(), ammonia = character()))
})

## For additional tests, refer to calculateFragments from PSMatch package
test_that("calculateFragments: Default behaviour without modifications", {
    
    ## Test 1: Default behavior without modifications
    sequence <- "PQR"
    result <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = NULL,
        variable_modifications = numeric(),
        max_mods = Inf,
        neutralLoss = defaultNeutralLoss(),
        verbose = FALSE
    )
    
    ## Check unique peptide without modifications
    expect_identical(unique(result$peptide), "PQR")
    })

test_that("calculateFragments: Behaviour with fixed modifications", {
    ## Test 2: Fixed modifications
    sequence <- "PQR"
    result <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = NULL,
        variable_modifications = NULL,
        max_mods = Inf,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE,
        modifications = NULL
    )
    fixed_modifications <- c(P = 79.966)
    result_fixed <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = fixed_modifications,
        variable_modifications = NULL,
        max_mods = 0,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE,
        modifications = NULL
    )
    
    ## Fixed modifications do not produce additional unique peptides
    expect_identical(unique(result_fixed$peptide), "PQR")
    expect_identical(nrow(result), nrow(result_fixed))
    
    ## Fixed modifications do change the fragment masses
    expect_false(all(result$mz == result_fixed$mz))
})

test_that("calculateFragments: Behaviour with variable modifications", {
    ## Test 3: Variable modifications only
    sequence <- "PQR"
    result <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = NULL,
        variable_modifications = NULL,
        max_mods = Inf,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE
    )
    
    fixed_modifications <- c(P = 79.966)
    result_fixed <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = fixed_modifications,
        variable_modifications = NULL,
        max_mods = 0,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE
    )
    
    variable_modifications <- c(P = 79.966, Q = 20, R = 10)
    max_mods <- 2
    result_var <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = NULL,
        variable_modifications = variable_modifications,
        max_mods = max_mods,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE
    )
    
    ## Calculate expected combinations
    expected_combinations <- 
        choose(3, 0) + choose(3, 1) + choose(3, 2) + choose(3, 3)
    
    ## Check if number of unique peptides matches expectations
    expect_equal(length(unique(result_var$peptide)), expected_combinations)
    
    ## Check if it's true in case there are less modifications than max_mods
    variable_modifications <- c(P = 79.966)
    max_mods <- 2
    result_var <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = NULL,
        variable_modifications = variable_modifications,
        max_mods = max_mods,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE
    )
    
    ## Calculate expected combinations
    expected_combinations <- choose(1, 0) + choose(1,1)
    
    ## Check if number of unique peptides matches expectations
    expect_equal(length(unique(result_var$peptide)), expected_combinations)
    
    ## Test 4: Fixed and variable modifications combined
    result_combined <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = fixed_modifications,
        variable_modifications = variable_modifications,
        max_mods = max_mods,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE
    )
    
    ## Check equal mass of variable mods fragments and no mods fragments
    expect_true(all(result$mz == result_var[result_var$peptide == "PQR","mz"]))
    
    ## Check equal mass of variable mods fragments and fixed mods fragments
    expect_true(all(result_fixed$mz == result_var[result_var$peptide == "[P]QR","mz"]))
})

test_that(".cumsumFragmentMasses: Behaviour with any modification", {
    ## Test4: Check behaviour of .cumsumFragmentMasses function
    
    ## Modifications used
    mods_forward <- c(P = 5, Q = 0, R = 7)
    mods_backward <- c(R = 7, Q = 0, P = 5)
    
    ## theoretical masses P = 15, Q = 25, R = 10)
    fragments_forward <- c(P = 15, Q = 40) ## representing cumsum forward ions
    fragments_backward <- c(R = 10, Q = 35) ## representing cumsum backward ions
    
    result_forward <- .cumsumFragmentMasses(mods_forward, fragments_forward)
    result_backward <- .cumsumFragmentMasses(mods_backward, fragments_backward)
    
    expect_identical(c(P = 20, Q = 45), result_forward)
    expect_identical(c(R = 17, Q = 42), result_backward)
})

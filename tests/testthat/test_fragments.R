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

    ## Deprecated parameters produce warnings
    expect_warning(
        calculateFragments("PQR", fixed_modifications = c(P = 2),
                           verbose = FALSE),
        "'fixed_modifications'.*deprecated"
    )
    expect_warning(
        calculateFragments("PQR", variable_modifications = c(P = 2),
                           verbose = FALSE),
        "'variable_modifications'.*deprecated"
    )

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
                 suppressWarnings(
                     calculateFragments(
                         "PQR",
                         fixed_modifications = c(C = 57.02146, Nterm = 229),
                         verbose = FALSE)),
                 check.attributes = FALSE, tolerance = 1e-5)

    ## neutral loss + nterm + cterm mod, rownames always differ
    tpqr$mz[3:7] <- tpqr$mz[3:7] - 100
    expect_equal(tpqr,
                 suppressWarnings(
                     calculateFragments(
                         "PQR",
                         fixed_modifications = c(C = 57.02146,
                                                 Nterm = 229,
                                                 Cterm = -100),
                         verbose = FALSE)),
                 check.attributes = FALSE, tolerance = 1e-5)

    expect_equal(ace,
                 calculateFragments("ACE", type = c("a", "b", "c", "x", "y", "z"),
                                    z = 2, neutralLoss = NULL,
                                    addCarbamidomethyl = FALSE,
                                    verbose = FALSE),
                 tolerance = 1e-5)
    expect_equal(ace[1:6,],
                 calculateFragments("ACE", type = letters[1:3], z = 2,
                                    addCarbamidomethyl = FALSE,
                                    verbose = FALSE),
                 tolerance = 1e-5)

    expect_error(calculateFragments("A"), "two or more residues")

    ## issue #200 (mz are not calculated correctly for terminal fixed_modifications
    ## and z > 1)
    p <- getAtomicMass()["p"]
    expect_equal(
        suppressWarnings(
            calculateFragments("AA", z = 2,
                               fixed_modifications = c(Nterm = 10),
                               type = "b"))$mz - p,
        (suppressWarnings(
            calculateFragments("AA", z = 1,
                               fixed_modifications = c(Nterm = 10),
                               type = "b"))$mz - p) / 2)
    expect_equal(
        suppressWarnings(
            calculateFragments("AA", z = 2, neutralLoss = NULL,
                               fixed_modifications = c(Cterm = 10),
                               type = "y"))$mz - p,
        (suppressWarnings(
            calculateFragments("AA", z = 1, neutralLoss = NULL,
                               fixed_modifications = c(Cterm = 10),
                               type = "y"))$mz - p) / 2)

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

test_that("calculateFragments: Default behaviour without modifications", {

    ## Default behavior without modifications
    sequence <- "PQR"
    result <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        addCarbamidomethyl = FALSE,
        neutralLoss = defaultNeutralLoss(),
        verbose = FALSE
    )

    ## Check unique peptide without modifications
    expect_identical(unique(result$peptide), "PQR")
})

test_that("calculateFragments: Default behaviour with positional modifications", {

    ## Default behavior with positional modifications in sequence string
    sequence <- "PQ[+10]R"
    result <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        addCarbamidomethyl = FALSE,
        neutralLoss = defaultNeutralLoss(),
        verbose = FALSE
    )

    ## Check unique peptide label preserved
    expect_identical(unique(result$peptide), "PQ[+10]R")
})

test_that("calculateFragments: Behaviour with fixed modifications", {
    ## Deprecated fixed_modifications path: still works, produces a warning
    sequence <- "PQR"
    result <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        addCarbamidomethyl = FALSE,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE
    )
    fixed_modifications <- c(P = 79.966)
    result_fixed <- suppressWarnings(
        calculateFragments(
            sequence = sequence,
            type = c("b", "y"),
            z = 1,
            fixed_modifications = fixed_modifications,
            neutralLoss = list(water = c(), ammonia = c()),
            verbose = FALSE
        )
    )

    ## Fixed modifications do not produce additional unique peptides
    expect_identical(unique(result_fixed$peptide), "PQR")
    expect_identical(nrow(result), nrow(result_fixed))

    ## Fixed modifications do change the fragment masses
    expect_false(all(result$mz == result_fixed$mz))
})

test_that("calculateFragments: Behaviour with variable modifications", {
    ## Deprecated variable_modifications path: still works, produces warnings
    sequence <- "PQR"
    result <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        addCarbamidomethyl = FALSE,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE
    )

    fixed_modifications <- c(P = 79.966)
    result_fixed <- suppressWarnings(
        calculateFragments(
            sequence = sequence,
            type = c("b", "y"),
            z = 1,
            fixed_modifications = fixed_modifications,
            neutralLoss = list(water = c(), ammonia = c()),
            verbose = FALSE
        )
    )

    variable_modifications <- c(P = 79.966, Q = 20, R = 10)
    max_mods <- 2
    result_var <- suppressWarnings(
        calculateFragments(
            sequence = sequence,
            type = c("b", "y"),
            z = 1,
            variable_modifications = variable_modifications,
            max_mods = max_mods,
            neutralLoss = list(water = c(), ammonia = c()),
            verbose = FALSE
        )
    )

    ## Calculate expected combinations
    expected_combinations <-
        choose(3, 0) + choose(3, 1) + choose(3, 2) + choose(3, 3)

    ## Check if number of unique peptides matches expectations
    expect_equal(length(unique(result_var$peptide)), expected_combinations)

    ## Check if it's true in case there are fewer modifications than max_mods
    variable_modifications <- c(P = 79.966)
    max_mods <- 2
    result_var <- suppressWarnings(
        calculateFragments(
            sequence = sequence,
            type = c("b", "y"),
            z = 1,
            variable_modifications = variable_modifications,
            max_mods = max_mods,
            neutralLoss = list(water = c(), ammonia = c()),
            verbose = FALSE
        )
    )

    ## Calculate expected combinations
    expected_combinations <- choose(1, 0) + choose(1, 1)

    ## Check if number of unique peptides matches expectations
    expect_equal(length(unique(result_var$peptide)), expected_combinations)

    ## Fixed and variable modifications combined
    result_combined <- suppressWarnings(
        calculateFragments(
            sequence = sequence,
            type = c("b", "y"),
            z = 1,
            fixed_modifications = fixed_modifications,
            variable_modifications = variable_modifications,
            max_mods = max_mods,
            neutralLoss = list(water = c(), ammonia = c()),
            verbose = FALSE
        )
    )

    ## Check equal mass of variable mods fragments and no mods fragments
    expect_true(all(result$mz == result_var[result_var$peptide == "PQR", "mz"]))

    ## Check equal mass of variable mods fragments and fixed mods fragments
    expect_true(
        all(result_fixed$mz == result_var[result_var$peptide == "[P]QR", "mz"])
    )
})

test_that("calculateFragments: addCarbamidomethyl parameter", {

    ## addCarbamidomethyl = FALSE: no extra mass on C-containing fragments
    result_no <- calculateFragments(
        "ACE",
        type = "b",
        neutralLoss = NULL,
        addCarbamidomethyl = FALSE,
        verbose = FALSE
    )

    ## addCarbamidomethyl = TRUE (default): adds 57.02146 to C-containing
    ## fragments
    result_cbm <- calculateFragments(
        "ACE",
        type = "b",
        neutralLoss = NULL,
        addCarbamidomethyl = TRUE,
        verbose = FALSE
    )

    ## b2 contains C (seq = "AC"), mass difference should be 57.02146
    b2_no  <- result_no[result_no$ion == "b2", "mz"]
    b2_cbm <- result_cbm[result_cbm$ion == "b2", "mz"]
    expect_equal(b2_cbm - b2_no, 57.02146, tolerance = 1e-5)

    ## b1 contains only A (seq = "A"), no difference expected
    b1_no  <- result_no[result_no$ion == "b1", "mz"]
    b1_cbm <- result_cbm[result_cbm$ion == "b1", "mz"]
    expect_equal(b1_cbm - b1_no, 0, tolerance = 1e-5)

    ## If [Carbamidomethyl] is already present in the sequence, it should not
    ## be added again
    result_present <- calculateFragments(
        "AC[Carbamidomethyl]E",
        type = "b",
        neutralLoss = NULL,
        addCarbamidomethyl = TRUE,
        verbose = FALSE
    )
    expect_equal(result_cbm$mz, result_present$mz, tolerance = 1e-5)

    ## Sequences without C are unaffected by addCarbamidomethyl
    result_pqr_false <- calculateFragments(
        "PQR", type = "b", neutralLoss = NULL,
        addCarbamidomethyl = FALSE, verbose = FALSE
    )
    result_pqr_true <- calculateFragments(
        "PQR", type = "b", neutralLoss = NULL,
        addCarbamidomethyl = TRUE, verbose = FALSE
    )
    expect_equal(result_pqr_false$mz, result_pqr_true$mz, tolerance = 1e-5)
})
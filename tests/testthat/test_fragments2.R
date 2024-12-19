## For additional tests, refer to calculateFragments from PSMatch package
test_that("calculateFragments2: Default behavior without modifications", {
    
    ## Test 1: Default behavior without modifications
    sequence <- "PQR"
    original_result <- calculateFragments(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        modifications = NULL,
        neutralLoss = defaultNeutralLoss(),
        verbose = FALSE
    )
    result <- calculateFragments2(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = NULL,
        variable_modifications = NULL,
        max_mods = Inf,
        neutralLoss = defaultNeutralLoss(),
        verbose = FALSE
    )
    
    ## Check unique peptide without modifications
    expect_equal(unique(result$peptide), "PQR")
    expect_equal(result[,-7], original_result)
    })

test_that("calculateFragments2: Behavior with fixed modifications", {
    ## Test 2: Fixed modifications
    sequence <- "PQR"
    result <- calculateFragments2(
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
    result_fixed <- calculateFragments2(
        sequence = sequence,
        type = c("b", "y"),
        z = 1,
        fixed_modifications = fixed_modifications,
        variable_modifications = NULL,
        max_mods = 0,
        neutralLoss = list(water = c(), ammonia = c()),
        verbose = FALSE
    )
    
    ## Fixed modifications do not produce additional unique peptides
    expect_equal(length(unique(result_fixed$peptide)), 1)
    
    ## Fixed modifications do change the fragment masses
    expect_false(all(result$mz == result_fixed$mz))
})

test_that("calculateFragments2: Behavior with variable modifications", {
    ## Test 3: Variable modifications only
    sequence <- "PQR"
    result <- calculateFragments2(
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
    result_fixed <- calculateFragments2(
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
    result_var <- calculateFragments2(
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
    expected_combinations <- choose(3, 0) + choose(3, 1) + choose(3, 2)
    
    ## Check if number of unique peptides matches expectations
    expect_equal(length(unique(result_var$peptide)), expected_combinations)
    
    ## Check if it's true in case there are less modifications than max_mods
    variable_modifications <- c(P = 79.966)
    max_mods <- 2
    result_var <- calculateFragments2(
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
    result_combined <- calculateFragments2(
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


setMethod("calculateFragments", c("character", "missing"),
          function(sequence, type = c("b", "y"), z = 1,
                   fixed_modifications = c(C = 57.02146),
                   variable_modifications = NULL,
                   max_mods = Inf,
                   neutralLoss = defaultNeutralLoss(),
                   verbose = TRUE) {
            l <- lapply(sequence, .calculateFragments,
                        type = type, z = z, modifications = modifications,
                        neutralLoss = neutralLoss, verbose = verbose)
            return(do.call(rbind, l))
        })
#' @title Calculate ions produced by fragmentation, including possible modifications
#' 
#' @param sequence character() providing a peptide sequence.
#' 
#' @param type character vector of target ions; possible values: 
#' c("a", "b", "c", "x", "y", "z"). Default is type = c("b", "y").
#' 
#' @param z numeric with a desired charge state; default is 1.
#' 
#' @param fixed_modifications A named numeric vector of used fixed modifications. 
#' The name must correspond to the one-letter-code of the modified amino acid 
#' and the numeric value must represent the mass that should be added to the 
#' original amino accid mass, default: Carbamidomethyl modifications = 
#' c(C = 57.02146). Use Nterm or Cterm as names for modifications that should 
#' be added to the amino respectively carboxyl-terminus.
#' 
#' @param variable_modifications A named numeric vector of variable modifications.
#' Depending on the maximum number of modifications (`max_mods`), all possible 
#' combinations are returned.
#' 
#' @param max_mods A numeric indicating the maximum number of variable modifications 
#' allowed on the sequence at once. Does not include fixed modifications. 
#' Default value is positive infinity.
#' 
#' @param neutralLoss list, it has to have two named elments, namely water and
#' ammonia that contain a character vector which type of neutral loss should be 
#' calculated. Currently neutral loss on the C terminal "Cterm", at the amino 
#' acids c("D", "E", "S", "T") for "water" (shown with an ⁠_⁠) an
#' d c("K", "N", "Q", "R") for "ammonia" (shown with an *) are supported.
#' 
#' @param verbose logical(1). If TRUE (default) the used modifications are printed.
#' 
#' @return A named list containing dataframes. The dataframes show all the 
#' ions produced by fragmentation with all possible combinations of modifications.
#' The used modifications are displayed in the element names of the list. Amino
#' acids within brackets thus include a variable modification. Fixed modifications
#' are not displayed in the names. 
#' 
#' @noRd
#' 
#' @examples
#' calculateFragments2(sequence = "ARGSHKATC", type = c("b", "y"), z = 1, 
#' fixed_modifications = c(C = 57), variable_modifications = c(S = 79, Y = 79, T = 79),
#' max_mods = 2)
#' 
calculateFragments2 <- function(sequence, 
                                type = c("b", "y"), 
                                z = 1,
                                fixed_modifications = c(C = 57.02146),
                                variable_modifications = NULL,
                                max_mods = Inf,
                                neutralLoss = defaultNeutralLoss(),
                                verbose = TRUE) {
    if (nchar(sequence) <= 1L) {
        stop("'sequence' has to have two or more residues.")
    }
    
    mod_combinations <- 
        .modificationPositions(sequence, 
                               variable_modifications, 
                               max_mods)
    type <- match.arg(type,
                      choices = c("a", "b", "c", "x", "y", "z"),
                      several.ok=TRUE)
    type <- sort(type)
    ## constants
    mass <- PSMatch::getAtomicMass()
    
    add <- c(a=-(mass["C"]+mass["O"]),            # + H - CO
             b=0,                                 # + H
             c=mass["N"]+3*mass["H"],             # + H + NH3
             x=mass["C"]+2*mass["O"],             # + CO + OH
             y=2*mass["H"]+mass["O"],             # + H2 + OH
             z=-(mass["N"]+mass["H"])+mass["O"])  # + NH + OH
    
    aa <- PSMatch::getAminoAcids()
    aamass <- setNames(aa$ResidueMass, aa$AA)
    
    ## replace default mass by masses with fixed modifications
    if (length(fixed_modifications)) {
        aamass[names(fixed_modifications)] <- 
            aamass[names(fixed_modifications)] + fixed_modifications
    }
    
    ## message used modifications
    if (verbose) {
        if (length(fixed_modifications)) {
            mods <-paste0(names(fixed_modifications),
                          "=",
                          fixed_modifications, 
                          collapse=", ")
        } else {
            mods <- "None"
        }
        if (length(variable_modifications)) {
            mods2 <- paste0(names(variable_modifications),
                            "=", 
                            variable_modifications, 
                            collapse=", ")
        } else {
            mods2 <- "None"
        }
        message("Fixed modifications used: ", mods, 
                "\nVariable modifications used: ", mods2)
    }
    
    ## split peptide sequence into aa
    fragment.seq <- strsplit(sequence, "")[[1]]
    fn <- length(fragment.seq)
    
    ## calculate cumulative mass starting at the amino-terminus (for a, b, c)
    amz <- cumsum(aamass[fragment.seq[-fn]])
    ## calculate cumulative mass starting at the carboxyl-terminus (for x, y, z)
    cmz <- cumsum(aamass[rev(fragment.seq[-1L])])
    
    ## calculate fragment mass (amino-terminus)
    tn <- length(amz)
    atype <- c("a", "b", "c") %in% type
    nat <- sum(atype)
    ## calculate fragment mass (carboxyl-terminus)
    ctype <- c("x", "y", "z") %in% type
    nct <- sum(ctype)
    
    ## devide by charge
    zn <- length(z)
    
    ## fragment seq (amino-terminus)
    aseq <- rep(rep(substring(sequence, rep(1L, fn - 1L),
                              1L:(fn - 1L)), each = zn), nat)
    
    ## fragment seq (carboxyl-terminus)
    cseq <- rep(rep(rev(substring(sequence, 2L:fn,
                                  rep(fn, fn - 1L))), each=zn), nct)
    
    ## add the variable modifications and apply steps above 
    amz_mod <- vector("list", length(mod_combinations))
    cmz_mod <- vector("list", length(mod_combinations))
    df <- vector("list", length(mod_combinations))
    
    for (i in 1:length(mod_combinations)) {
        amz_mod[[i]] <- .cumsumFragmentMasses(mod_combinations[[i]], amz)
        cmz_mod[[i]] <- .cumsumFragmentMasses(rev(mod_combinations[[i]]), cmz)
        
        amz_mod[[i]] <- rep(amz_mod[[i]], nat) + rep(add[1:3][atype], each=tn)
        cmz_mod[[i]] <- rep(cmz_mod[[i]], nct) + rep(add[4:6][ctype], each=tn)
        
        amz_mod[[i]] <- rep(amz_mod[[i]], each = zn)/z
        cmz_mod[[i]] <- rep(cmz_mod[[i]], each = zn)/z
        
        ## add protons (H+)
        amz_mod[[i]] <- amz_mod[[i]] + mass["p"]
        cmz_mod[[i]] <- cmz_mod[[i]] + mass["p"]
    }
    
    ## fragment str (amino-terminus)
    atype <- rep(c("a", "b", "c")[atype], each = tn * zn)
    pos <- rep(1L:tn, each = zn)
    if (length(atype)) {
        aion <- paste0(atype, pos)
    } else {
        aion <- character()
    }
    
    ## fragment str (carboxyl-terminus)
    ctype <- rep(c("x", "y", "z")[ctype], each = tn * zn)
    if (length(ctype)) {
        cion <- paste0(ctype, pos)
    } else {
        cion <- character()
    }
    
    ## generate unique dataframe with all fragments and modifications 
    for (i in 1:length(mod_combinations)) {
        df[[i]] <- data.frame(mz = c(amz_mod[[i]], cmz_mod[[i]]),
                              ion = c(aion, cion),
                              type = c(atype, ctype),
                              pos = pos,
                              z = z,
                              seq = c(aseq, cseq),
                              stringsAsFactors = FALSE)
        df[[i]] <- .neutralLoss(df[[i]],
                                water = neutralLoss$water,
                                ammonia = neutralLoss$ammonia)
        df[[i]] <- .terminalModifications(df[[i]],
                                          modifications = fixed_modifications)
        rownames(df[[i]]) <- NULL
        non_zero <- mod_combinations[[i]] != 0
        names(mod_combinations[[i]])[non_zero] <- 
            paste0("[", (names(mod_combinations[[i]])[non_zero]), "]")
        df[[i]][["peptide"]] <- paste(names(mod_combinations[[i]]),
                                      collapse = "")
    }
    
    df <- do.call(rbind, df)
    rownames(df) <- NULL
    df
}

#' @title Generates list of possible combinations of modifications
#' 
#' @param sequence Character. A peptide sequence that may have modifications or not 
#' 
#' @param fixed_modifications Named numeric. Specifies which fixed modifications are used
#' 
#' @param variable_modifications Named numeric. Specifies which variable modifications are used
#' 
#' @param max_mods Numeric. Indicates how many modifications can be applied at once.
#' 
#' @return list with all possible combinations of modifications
#' 
#' @noRd
#' 
#' @examples
#' .modificationPositions("ARGHKA", variable_modifications = c(A = 4, K = 5, S = 8), max_mods = 3)
#' 
.modificationPositions <- function(sequence,
                                   variable_modifications = NULL,
                                   max_mods = Inf) {
    
    sequence_split <- strsplit(sequence, "")[[1]]
    
    if (length(variable_modifications)) {
        modifiable_positions_var <-
            which(sequence_split %in% names(variable_modifications))
        l <- length(modifiable_positions_var)
    } else {l <- NULL}
    
    ## take the maximum amount of modifications possible 
    max_mods <- min(max_mods, l)
    
    if (max_mods > 0 & length(l)) {
        mod_combinations <- lapply(0:max_mods, function(n) {
            if (n == 0) {
                list(integer(0))
            } else if (l == 1) {
                list(modifiable_positions_var)
            } else {
                combn(modifiable_positions_var, n, simplify = FALSE)
            }
        })
        
        mod_combinations <- do.call(c, mod_combinations)
        lapply(mod_combinations, function(comb) {
            mods <- rep(0, nchar(sequence))
            mods[comb] <- unname(variable_modifications[sequence_split[comb]])
            names(mods) <- sequence_split
            mods
        })
    } else {
        mods <- rep(0, nchar(sequence))
        names(mods) <- sequence_split
        list(mods)
    }
}

# .modificationPositions <- function(sequence,
#                                    variable_modifications = NULL,
#                                    max_mods = 2) {
#     
#     sequence_split <- strsplit(sequence, "")[[1]]
#     
#     if (length(variable_modifications)) {
#         modifiable_positions_var <-
#             which(sequence_split %in% names(variable_modifications))
#         l <- length(modifiable_positions_var)
#     } else {l <- NULL}
#     
#     ## take the maximum amount of modifications possible 
#     max_mods <- min(max_mods, l)
#     
#     if (max_mods > 0 & length(l)) {
#         mod_combinations <- lapply(0:max_mods, function(n) {
#             if (n == 0) {
#                 list(integer(0))
#             } else if (l == 1) {
#                 list(modifiable_positions_var)
#             } else {
#                 combn(modifiable_positions_var, n, simplify = FALSE)
#             }
#         })
#         
#         mod_combinations <- do.call(c, mod_combinations)
#         lapply(mod_combinations, function(comb) {
#             mods <- rep(0, nchar(sequence))
#             mods[comb] <- unname(variable_modifications[sequence_split[comb]])
#             names(mods) <- sequence_split
#             mods
#         })
#     } else {
#         mods <- rep(0, nchar(sequence))
#         names(mods) <- sequence_split
#         mods
#     }
# }
#' @title Sums the modifications cumulatively to the sequence 
#' 
#' @param modificationCombination Named numeric. Element of list of `.modification_positions`
#' 
#' @param fragmentMasses Named numeric. Fragment masses for a given sequence
#' 
#' @return A named numeric indicating the fragment masses for a given sequence
#' with the assigned combination of modifications
#' 
#' @noRd
#' 
#' @examples
#' first <- c(A = 4, R = 0, G = 0, H = 0, K = 0, A = 6)
#' second <- c(A = 75, R = 231, G = 288, H = 425, K = 553, A = 580)
#' .cumsumFragmentMasses(first, second)

.cumsumFragmentMasses <- function(modificationCombination, fragmentMasses) {
    
    expanded_modificationCombination <- rep(0, length(fragmentMasses))
    cumulative_sum <- 0
    
    for (i in seq_along(fragmentMasses)) {
        cumulative_sum <- cumulative_sum + modificationCombination[i]
        expanded_modificationCombination[i] <- cumulative_sum
    }
    
    fragmentMasses + expanded_modificationCombination
}

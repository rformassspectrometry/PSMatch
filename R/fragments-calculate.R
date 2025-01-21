#' @title Calculate ions produced by fragmentation with variable modifications
#' 
#' @aliases calculateFragments modificationPositions defaultNeutralLoss calculateFragments,character,missing-method
#'
#' @name calculateFragments
#' 
#' @description
#'
#' This method calculates a-, b-, c-, x-, y- and z-ions produced by
#' fragmentation.
#'
#' Available methods
#'
#' - The default method with signature `sequence = "character"` and
#'   `object = "missing"` calculates the theoretical fragments for a
#'   peptide sequence. It returns a `data.frame` with the columns
#'   `mz`, `ion`, `type`, `pos`, `z`, `seq` and `peptide`.
#'
#' - Additional method can be defined that will adapt their behaviour
#'   based on spectra defined in `object`. See for example the MSnbase
#'   package that implements a method for objects of class
#'   `Spectrum2`.
#'   
#' @param sequence `character()` providing a peptide sequence.
#' 
#' @param type `character` vector of target ions; possible values: 
#' `c("a", "b", "c", "x", "y", "z")`. Default is `type = c("b", "y")`.
#' 
#' @param z numeric with a desired charge state; default is 1.
#' 
#' @param fixed_modifications A named `numeric` vector of used fixed modifications. 
#' The name must correspond to the one-letter-code of the modified amino acid 
#' and the numeric value must represent the mass that should be added to the 
#' original amino accid mass, default: Carbamidomethyl modifications = 
#' c(C = 57.02146). Use Nterm or Cterm as names for modifications that should 
#' be added to the amino respectively carboxyl-terminus.
#' 
#' @param variable_modifications A named `numeric` vector of variable modifications.
#' Depending on the maximum number of modifications (`max_mods`), all possible 
#' combinations are returned.
#' 
#' @param max_mods A numeric indicating the maximum number of variable modifications 
#' allowed on the sequence at once. Does not include fixed modifications. 
#' Default value is positive infinity.
#' 
#' @param neutralLoss `list`, it has to have two named elments,
#'     namely `water` and `ammonia` that contain a `character` vector
#'     which type of neutral loss should be calculated.  Currently
#'     neutral loss on the C terminal `"Cterm"`, at the amino acids
#'     `c("D", "E", "S", "T")` for `"water"` (shown with an `_`) and
#'     `c("K", "N", "Q", "R")` for `"ammonia"` (shown with an `*`) are
#'     supported.
#'
#'     There is a helper function `defaultNeutralLoss()` that returns
#'     the correct list. It has two arguments `disableWaterLoss` and
#'     `disableAmmoniaLoss` to remove single neutral loss options. See
#'     the example section for use cases.
#' 
#' @param verbose `logical(1)`. If `TRUE` (default) the used modifications are printed.
#' 
#' @param modifications Named `numeric()`. Deprecated modifications parameter.
#' Will override `fixed_modifications` but is set to `NULL` by default. Please 
#' refrain from using it, opt for `fixed_modifications` instead.
#' 
#' @return A `data.frame` showing all the 
#' ions produced by fragmentation with all possible combinations of modifications.
#' The used variable modifications are displayed in the `peptide` column through the
#' use of amino acids followed by the modification within brackets. 
#' Fixed modifications are not displayed. 
#' 
#' @author Sebastian Gibb <mail@sebastiangibb.de>
#' 
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @importFrom ProtGenerics calculateFragments
#'
#' @exportMethod calculateFragments
#' 
#' @examples
#' ## General use
#' calculateFragments(sequence = "ARGSHKATC", type = c("b", "y"), z = 1, 
#' fixed_modifications = c(C = 57), variable_modifications = c(S = 79, Y = 79, T = 79),
#' max_mods = 2)
#' 
#' ## calculate fragments for ACE with default modification
#' calculateFragments("ACE", fixed_modifications = c(C = 57.02146))
#'
#' #' ## calculate fragments for ACE with an added variable modification
#' calculateFragments("ACE", variable_modifications = c(A = 43.25))
#' 
#' ## calculate fragments for ACE with an added N-terminal modification
#' calculateFragments("ACE", fixed_modifications = c(C = 57.02146, Nterm = 229.1629))
#'
#' ## calculate fragments for ACE without any modifications
#' calculateFragments("ACE", fixed_modifications = NULL)
#'
#' calculateFragments("VESITARHGEVLQLRPK",
#'                    type = c("a", "b", "c", "x", "y", "z"),
#'                    z = 1:2)
#'
#' ## neutral loss
#' defaultNeutralLoss()
#'
#' ## disable water loss on the C terminal
#' defaultNeutralLoss(disableWaterLoss="Cterm")
#'
#' ## real example
#' calculateFragments("PQR")
#' calculateFragments("PQR",
#'                    neutralLoss=defaultNeutralLoss(disableWaterLoss="Cterm"))
#' calculateFragments("PQR",
#'                    neutralLoss=defaultNeutralLoss(disableAmmoniaLoss="Q"))
#'
#' ## disable neutral loss completely
#' calculateFragments("PQR", neutralLoss=NULL)
#'  
setMethod("calculateFragments", c("character", "missing"),
          function(sequence, type = c("b", "y"), z = 1,
                   fixed_modifications = c(C = 57.02146),
                   variable_modifications = numeric(),
                   max_mods = Inf,
                   neutralLoss = defaultNeutralLoss(),
                   verbose = TRUE, 
                   modifications = NULL) {
              l <- lapply(sequence, .calculateFragments,
                          type = type, z = z, 
                          fixed_modifications = fixed_modifications,
                          variable_modifications = variable_modifications,
                          max_mods = Inf,
                          neutralLoss = neutralLoss, 
                          verbose = verbose,
                          modifications = modifications)
              return(do.call(rbind, l))
          })


#' calculate fragments from a peptide sequence
#'
#' @param sequence character vector of length 1
#'
#' @param type could be c("a", "b", "c", "x", "y", "z")
#'
#' @param z charge
#'
#' @param fixed_modifications A named `numeric` vector of used fixed modifications. 
#' The name must correspond to the one-letter-code of the modified amino acid 
#' and the numeric value must represent the mass that should be added to the 
#' original amino accid mass, default: Carbamidomethyl modifications = 
#' c(C = 57.02146). Use Nterm or Cterm as names for modifications that should 
#' be added to the amino respectively carboxyl-terminus.
#' 
#' @param variable_modifications A named `numeric` vector of variable modifications.
#' Depending on the maximum number of modifications (`max_mods`), all possible 
#' combinations are returned.
#' 
#' @param max_mods A numeric indicating the maximum number of variable modifications 
#' allowed on the sequence at once. Does not include fixed modifications. 
#' Default value is positive infinity.
#' 
#' @param neutralLoss `list`, it has to have two named elments,
#'     namely `water` and `ammonia` that contain a `character` vector
#'     which type of neutral loss should be calculated.  Currently
#'     neutral loss on the C terminal `"Cterm"`, at the amino acids
#'     `c("D", "E", "S", "T")` for `"water"` (shown with an `_`) and
#'     `c("K", "N", "Q", "R")` for `"ammonia"` (shown with an `*`) are
#'     supported.
#'
#'     There is a helper function `defaultNeutralLoss()` that returns
#'     the correct list. It has two arguments `disableWaterLoss` and
#'     `disableAmmoniaLoss` to remove single neutral loss options. See
#'     the example section for use cases.
#' 
#' @param verbose `logical(1)`. If `TRUE` (default) the used modifications are printed.
#' 
#' @param modifications Named `numeric()`. Deprecated modifications parameter.
#' Will override `fixed_modifications` but is set to `NULL` by default. Please 
#' refrain from using it, opt for `fixed_modifications` instead.
#' 
#' @importFrom stats setNames
#'
#' @noRd
.calculateFragments <- function(sequence, 
                                type = c("b", "y"), 
                                z = 1,
                                fixed_modifications = c(C = 57.02146),
                                variable_modifications = numeric(),
                                max_mods = Inf,
                                neutralLoss = defaultNeutralLoss(),
                                verbose = TRUE,
                                modifications = NULL) {
    if (nchar(sequence) <= 1L) {
        stop("'sequence' has to have two or more residues.")
    }
    
    if (!is.null(modifications)) {
        warning("'modifications' is deprecated, please use 'fixed_modifications' instead.")
        fixed_modifications <- modifications
    }
    
    ## split peptide sequence into aa
    fragment.seq <- strsplit(sequence, "")[[1]]
    fn <- length(fragment.seq)
    
    mod_combinations <- 
        .modificationPositions(fragment.seq, 
                               variable_modifications, 
                               max_mods)
    type <- match.arg(type,
                      choices = c("a", "b", "c", "x", "y", "z"),
                      several.ok=TRUE)
    type <- sort(type)
    ## constants
    mass <- getAtomicMass()
    ## according to Table 1 of:
    ## Johnson, R. S., Martin, S. A., Biemann, K., Stults, J. T., and
    ## Watson, J. T. (1987).
    ## Novel fragmentation process of peptides by collision-induced
    ## decomposition in a tandem mass spectrometer: differentiation of leucine
    ## and isoleucine.
    ## Analytical Chemistry, 59(21), 2621-2625.
    ## https://doi.org/10.1021/ac00148a019
    ##
    ## a proton (H+) is added later
    ## (after calculation of the different charge states)
    add <- c(a=-(mass["C"]+mass["O"]),            # + H - CO
             b=0,                                 # + H
             c=mass["N"]+3*mass["H"],             # + H + NH3
             x=mass["C"]+2*mass["O"],             # + CO + OH
             y=2*mass["H"]+mass["O"],             # + H2 + OH
             z=-(mass["N"]+mass["H"])+mass["O"])  # + NH + OH
    
    aa <- getAminoAcids()
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
            paste0(names(mod_combinations[[i]])[non_zero],
                   "[",
                   mod_combinations[[i]][non_zero],
                   "]")
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
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#' 
#' @importFrom utils combn
#'
#' @noRd
#' 
#' @examples
#' .modificationPositions("ARGHKA", variable_modifications = c(A = 4, K = 5, S = 8), max_mods = 3)
#' 
.modificationPositions <- function(fragment.seq,
                                   variable_modifications = numeric(),
                                   max_mods = Inf) {
    modifiable_positions_var <- 
        which(fragment.seq %in% names(variable_modifications))

    l <- length(modifiable_positions_var)

    ## take the maximum amount of modifications possible
    max_mods <- min(max_mods, l)

    if (!length(variable_modifications) || max_mods <= 0)
        return(
            list(setNames(integer(length(fragment.seq)), fragment.seq))
        )

    .mod <- function(cmb,
                     seq_split = fragment.seq,
                     var_mods = variable_modifications) {
        m <- setNames(integer(length(seq_split)), seq_split)
        m[cmb] <- var_mods[seq_split[cmb]]
        m
    }

    c(
        list(setNames(integer(length(fragment.seq)), fragment.seq)),
        if (length(modifiable_positions_var) == 1)
            lapply(modifiable_positions_var, .mod)
        else
            unlist(
                lapply(seq_len(max_mods),
                    function(n)combn(
                        modifiable_positions_var, n,
                        FUN = .mod,
                        simplify = FALSE
                    )
                ),
                recursive = FALSE
            )
    )
}

.cumsumFragmentMasses <- function(modificationCombination, fragmentMasses) {
    
    modificationCombination <- 
        modificationCombination[-NROW(modificationCombination)]
    
    fragmentMasses + cumsum(modificationCombination)
}

#' adds neutral loss to data.frame generated by .calculateFragments
#' 
#' @param df data.frame generated by. calculateFragments
#' 
#' @return data.frame neutral loss rows added
#' 
#' @noRd
.neutralLoss <- function(df,
                         water = c("Cterm", "D", "E", "S", "T"),
                         ammonia = c("K", "N", "Q", "R")) {
    ## see "Low energy peptide fragmentation pathways" by Hugh-G. Patterton, Ph.D.
    ## http://cbio.ufs.ac.za/fgap/download/fragmentation_review.pdf
    ## see also discussion #47: https://github.com/lgatto/MSnbase/issues/47
    
    ## constants
    mass <- getAtomicMass()
    
    widx <- double()
    aidx <- double()
    
    .removeNeutralLoss <- function(df, idx, mass, ion) {
        if (length(idx)) {
            loss <- df[idx, ]
            loss[, c("ion", "type")] <- paste0(c(loss$ion, loss$type), ion)
            loss$mz <- loss$mz - mass / loss$z
            rbind(df, loss)
        } else {
            df
        }
    }
    
    if (length(water)) {
        ## N-term D/E, internal S/T
        rules <- c(D = "^D.", E = "^E.", S = ".S.", T = ".T.")
        rules <- rules[intersect(c("D", "E", "S", "T"), water)]
        
        if (length(rules)) {
            widx <- grep(paste0(rules, collapse = "|"), df$seq)
        }
        
        ## C-term COOH (all x, y, z fragments)
        if ("Cterm" %in% water) {
            widx <- unique(c(widx, grep("[xyz]", df$type)))
        }
    }
    
    if (length(ammonia)) {
        ## N-term/internal K/N/Q, internal R
        rules <- c(K = "^.*K.", N = "^.*N.", Q = "^.*Q.", R = ".R.")
        rules <- rules[intersect(c("K", "N", "Q", "R"), ammonia)]
        
        if (length(rules)) {
            aidx <- grep(paste0(rules, collapse="|"), df$seq)
        }
    }
    
    if (length(widx)) {
        df <- .removeNeutralLoss(df, idx = widx, mass = 2*mass["H"]+mass["O"], ion = "_")
    }
    if (length(aidx)) {
        df <- .removeNeutralLoss(df, idx = aidx, mass = mass["N"]+3*mass["H"], ion = "*")
    }
    df
}

#' adds nterm/cterm modifications to data.frame generated by
#' .calculateFragments should be used after .neutralLoss
#'
#' @param df data.frame generated by. calculateFragments
#'
#' @return modified data.frame
#'
#' @noRd
.terminalModifications <- function(df, modifications) {
    
    if ("Nterm" %in% names(modifications)) {
        isABC <- grep("[abc]", df$type)
        
        if (length(isABC)) {
            df$mz[isABC] <- df$mz[isABC] + modifications["Nterm"] / df$z[isABC]
        }
    }
    
    if ("Cterm" %in% names(modifications)) {
        isXYZ <- grep("[xyz]", df$type)
        
        if (length(isXYZ)) {
            df$mz[isXYZ] <- df$mz[isXYZ] + modifications["Cterm"] / df$z[isXYZ]
        }
    }
    
    df
}

#' default neutral loss argument for calculateFragments
#'
#' @param disableWaterLoss character, which loss should not calculated
#'
#' @param disableAmmoniaLoss character, which loss should not
#'     calculated
#'
#' @export
#'
#' @noRd
defaultNeutralLoss <- function(disableWaterLoss = NULL, disableAmmoniaLoss = NULL) {
    list(water = setdiff(c("Cterm", "D", "E", "S", "T"), disableWaterLoss),
         ammonia = setdiff(c("K", "N", "Q", "R"), disableAmmoniaLoss))
}
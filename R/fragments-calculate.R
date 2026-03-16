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
#' @param sequence `character()` providing a peptide sequence. If positional
#' modifications are included in the sequence, variable modifications may not be
#' used. See examples below for more detail.
#'
#' @param type `character` vector of target ions; possible values:
#' `c("a", "b", "c", "x", "y", "z")`. Default is `type = c("b", "y")`.
#'
#' @param z numeric with a desired charge state; default is 1.
#'
#' @param fixed_modifications Deprecated parameter. Please use
#' [PTMods::addFixedModifications()] to generate sequences with positional
#' modifications instead. Named `numeric` or `character`. If a `character` is
#' given, values must be in UniMod name or UniMod ID format
#' (e.g. `"Phospho"`, `"UNIMOD:21"`). The annotation style of
#' the values is preserved in the output. Specifies which fixed
#' modifications are applied to which amino acids.
#'
#' @param variable_modifications Deprecated parameter. Please use
#' [PTMods::addVariableModifications()] to generate sequences with positional
#' modifications instead. Named `numeric` or `character`. If a `character`
#' is given, values must be in UniMod name or UniMod ID format
#' (e.g. `"Phospho"`, `"UNIMOD:21"`). The annotation style of
#' the values is preserved in the output. Specifies which
#' variable modifications are used on which amino acids.
#'
#' @param max_mods Deprecated parameter. Please use in combination with
#' PTMods::addVariableModifications() instead. A numeric indicating the
#' maximum number of variable modifications
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
#' @param verbose `logical(1)`. Deprecated parameter.
#' If `TRUE` (default) the used modifications are printed.
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
#'
#' ## No modifications
#' calculateFragments("ACE")
#'
#' ## Multiple ion types and charge states
#' calculateFragments("ACE",
#'                    type = c("a", "b", "c", "x", "y", "z"),
#'                    z = 1:2)
#'
#' ## Positional modification written directly in the sequence string
#' ## The annotation style must be supported by PTMods::convertAnnotation
#' calculateFragments("A[+43.25]CE")
#' calculateFragments("T[Phospho]CE")
#' calculateFragments("T[UNIMOD:21]C[Carbamidomethyl]E")
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
#' ## Recommended workflow: use PTMods functions to produce positional sequences
#' ## before calling calculateFragments.
#'
#' ## Fixed modification (Carbamidomethyl on C) using addFixedModifications
#' seq_fixed <- PTMods::addFixedModifications("ACE",
#'                                            fixedModifications = c(C = 57.02))
#' calculateFragments(seq_fixed)
#'
#' ## Fixed modification including N-terminus using addFixedModifications
#' seq_nterm <- PTMods::addFixedModifications(
#'     "ACE",
#'     fixedModifications = c(C = 57.02, Nterm = 229.16))
#' calculateFragments(seq_nterm)
#'
#' ## Variable modification (delta mass on A) using addVariableModifications
#' seq_var <- PTMods::addVariableModifications("ACE",
#'                                             variableModifications = c(A = 43.25))
#' calculateFragments(seq_var)
#'
#' ## Both fixed and variable modifications using addModifications
#' seq_mods <- PTMods::addModifications("ARGSHKATC",
#'                                      fixedModifications = c(C = 57),
#'                                      variableModifications = c(S = 79, T = 79),
#'                                      maxMods = 2)
#' calculateFragments(seq_mods)
setMethod("calculateFragments", c("character", "missing"),
          function(sequence, type = c("b", "y"), z = 1,
                   fixed_modifications = NULL,
                   variable_modifications = NULL,
                   max_mods = Inf,
                   neutralLoss = defaultNeutralLoss(),
                   verbose = TRUE,
                   modifications = NULL) {

        if (!is.null(modifications)) {
            warning("'modifications' is deprecated,
                please use 'PTMods::addFixedModifications' instead.")
            fixed_modifications <- modifications
        }

        if (!is.null(fixed_modifications)) {
            warning("'fixed_modifications' is deprecated,
                please use 'PTMods::addFixedModifications' instead.")
            fixed_modifications <- modifications
        }

        if (!is.null(variable_modifications)) {
            warning("'variable_modifications' is deprecated,
                please use 'PTMods::addVariableModifications' instead.")
            fixed_modifications <- modifications
        }

        if (!is.null(fixed_modifications) | !is.null(variable_modifications)) {
            sequence <- sequence |>
                PTMods::addFixedModifications(
                    fixedModifications = fixed_modifications) |>
                PTMods::addVariableModifications(
                    variable_modifications, maxMods = max_mods)
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
                    collapse = ", ")
            } else {
                mods2 <- "None"
            }

            if (length(fixed_modifications) | length(variable_modifications)) {
                message("Fixed modifications used: ", mods,
                    "\nVariable modifications used: ", mods2)
            }
        }

        l <- lapply(sequence, .calculateFragments,
            type = type, z = z, neutralLoss = neutralLoss)

        return(do.call(rbind, l))
    }
)

#' calculate fragments from a peptide sequence
#'
#' @param sequence character vector of length 1
#'
#' @param type could be c("a", "b", "c", "x", "y", "z")
#'
#' @param z charge
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
#' @importFrom stats setNames
#'
#' @noRd
.calculateFragments <- function(sequence,
                                type = c("b", "y"),
                                z = 1,
                                neutralLoss = defaultNeutralLoss()) {

    initial_sequence <- sequence
    sequence <- convertAnnotation(sequence)

    parsed_modifications <- PTMods:::.parseModifiedSequence(sequence)
    canonical_sequence <- PTMods::getCanonicalSequence(sequence)

    if (nchar(canonical_sequence) <= 1L) {
        stop("'sequence' has to have two or more residues.")
    }

    ## split peptide sequence into aa
    fragment.seq <- strsplit(canonical_sequence, "")[[1]]
    fn <- length(fragment.seq)

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

    ## calculate cumulative mass starting at the amino-terminus (for a, b, c)
    amz <- cumsum(parsed_modifications[-fn]) + cumsum(aamass[fragment.seq[-fn]])
    ## calculate cumulative mass starting at the carboxyl-terminus (for x, y, z)
    cmz <- cumsum(rev(parsed_modifications[-1L])) + cumsum(aamass[rev(fragment.seq[-1L])])

    ## calculate fragment mass (amino-terminus)
    tn <- length(amz)
    atype <- c("a", "b", "c") %in% type
    nat <- sum(atype)
    amz <- rep(amz, nat) + rep(add[1:3][atype], each=tn)

    ## calculate fragment mass (carboxyl-terminus)
    ctype <- c("x", "y", "z") %in% type
    nct <- sum(ctype)
    cmz <- rep(cmz, nct) + rep(add[4:6][ctype], each=tn)

    ## devide by charge
    zn <- length(z)
    amz <- rep(amz, each = zn)/z
    cmz <- rep(cmz, each = zn)/z

    ## add protons (H+)
    amz <- amz + mass["p"]
    cmz <- cmz + mass["p"]

    ## fragment seq (amino-terminus)
    aseq <- rep(rep(substring(canonical_sequence, rep(1L, fn - 1L),
                              1L:(fn - 1L)), each = zn), nat)

    ## fragment seq (carboxyl-terminus)
    cseq <- rep(rep(rev(substring(canonical_sequence, 2L:fn,
                                  rep(fn, fn - 1L))), each=zn), nct)

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

    df <- data.frame(mz = c(amz, cmz),
                      ion = c(aion, cion),
                      type = c(atype, ctype),
                      pos = pos,
                      z = z,
                      seq = c(aseq, cseq),
                      stringsAsFactors = FALSE)
    df <- .neutralLoss(df,
        water = neutralLoss$water,
        ammonia = neutralLoss$ammonia)
    df <- .addTerminalMasses(df,
        parsedModifications = parsed_modifications)

    df$peptide <- initial_sequence
    rownames(df) <- NULL
    df
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
#' @return termini-modified data.frame
#'
#' @noRd
.addTerminalMasses <- function(df, parsedModifications) {

    nterm <- attr(parsedModifications, "Nterm")
    cterm <- attr(parsedModifications, "Cterm")

    if (!is.na(nterm)) {
        isABC <- grep("[abc]", df$type)

        if (length(isABC)) {
            df$mz[isABC] <- df$mz[isABC] + nterm / df$z[isABC]
        }
    }

    if (!is.na(cterm)) {
            isXYZ <- grep("[xyz]", df$type)

            if (length(isXYZ)) {
                df$mz[isXYZ] <- df$mz[isXYZ] + cterm / df$z[isXYZ]
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
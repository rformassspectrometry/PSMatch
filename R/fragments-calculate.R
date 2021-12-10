#' @title Calculate ions produced by fragmentation
#'
#' @aliases defaultNeutralLoss calculateFragments calculateFragments,character,missing-method
#'
#' @name calculateFragments
#' 
#' @description
#'
#' These method calculates a-, b-, c-, x-, y- and z-ions produced by
#' fragmentation.
#'
#' Available methods
#' 
#' - The default method with signature `sequence = "character"` and
#'   `object = "missing"` calculates the theoretical fragments for a
#'   peptide sequence. It returns a `data.frame` with the columns
#'   `mz`, `ion`, `type`, `pos`, `z` and `seq`.
#'
#' - Additional method can be defined that will adapt their behaviour
#'   based on spectra defined in `object`. See for example the MSnbase
#'   package that implements a method for objects of class
#'   `Spectrum2`.
#'
#' @return The methods with `oject = "missing"` returns a
#'     `data.frame`.
#' 
#' @param sequence character() providing a peptide sequence.
#' 
#' @param type `character` vector of target ions; possible values:
#'     `c("a", "b", "c", "x", "y", "z")`. Default is `type = c("b",
#'     "y")`.
#' 
#' @param z `numeric` with desired charge state; default is 1.
#' 
#' @param modifications A named `numeric` vector of used
#'     modifications. The name must correspond to the one-letter-code
#'     of the modified amino acid and the `numeric` value must
#'     represent the mass that should be added to the original amino
#'     accid mass, default: Carbamidomethyl `modifications = c(C =
#'     57.02146)`. Use `Nterm` or `Cterm` as names for modifications
#'     that should be added to the amino respectively
#'     carboxyl-terminus.
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
#' @param verbose `logical(1)`. If `TRUE` (default) the used
#'     modifications are printed.
#'
#' @author Sebastian Gibb <mail@sebastiangibb.de>
#' 
#' @importFrom ProtGenerics calculateFragments
#'
#' @exportMethod calculateFragments
#'
#' @examples
#'
#' ## calculate fragments for ACE with default modification
#' calculateFragments("ACE", modifications = c(C = 57.02146))
#'
#' ## calculate fragments for ACE with an addition N-terminal modification
#' calculateFragments("ACE", modifications = c(C = 57.02146, Nterm = 229.1629))
#'
#' ## calculate fragments for ACE without any modifications
#' calculateFragments("ACE", modifications = NULL)
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
setMethod("calculateFragments", c("character", "missing"),
          function(sequence, type = c("b", "y"), z = 1,
                   modifications = c(C = 57.02146),
                   neutralLoss = defaultNeutralLoss(),
                   verbose = TRUE) {
            l <- lapply(sequence, .calculateFragments,
                        type = type, z = z, modifications = modifications,
                        neutralLoss = neutralLoss, verbose = verbose)
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
#' @param modifications a named (amino acid one-letter-code; upper
#'     case) vector of mass that should be added (default:
#'     Carbamidomethyl 57.02146 is added to Cystein (C) and results in
#'     160.030649).
#' 
#' @param neutralLoss list, currently water and ammonia loss are supported
#' 
#' @param verbose verbose output?
#'
#' @importFrom stats setNames
#' 
#' @noRd
.calculateFragments <- function(sequence, type = c("b", "y"), z = 1,
                                modifications = c(C = 57.02146),
                                neutralLoss = defaultNeutralLoss(),
                                verbose = TRUE) {
    if (nchar(sequence) <= 1L) {
        stop("'sequence' has to have two or more residues.")
    }

    type <- match.arg(type, choices=c("a", "b", "c", "x", "y", "z"), several.ok=TRUE)
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

    ## replace default mass by modifications
    if (length(modifications)) {
        aamass[names(modifications)] <- aamass[names(modifications)] + modifications
    }

    if (verbose) {
        if (length(modifications)) {
            mods <- paste0(names(modifications), "=", modifications, collapse=", ")
        } else {
            mods <- "None"
        }
        message("Modifications used: ", mods)
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
    aseq <- rep(rep(substring(sequence, rep(1L, fn - 1L),
                              1L:(fn - 1L)), each = zn), nat)

    ## fragment seq (carboxyl-terminus)
    cseq <- rep(rep(rev(substring(sequence, 2L:fn,
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

    df <- data.frame(mz = c(amz, cmz),
                     ion = c(aion, cion),
                     type = c(atype, ctype),
                     pos = pos,
                     z = z,
                     seq = c(aseq, cseq),
                     stringsAsFactors = FALSE)
    df <- .neutralLoss(df, water = neutralLoss$water, ammonia = neutralLoss$ammonia)
    df <- .terminalModifications(df, modifications = modifications)
    rownames(df) <- NULL
    df
}

#' adds neutral loss to data.frame generated by .calculateFragments
#' @param df data.frame generated by. calculateFragments
#' @return data.frame neutral loss rows added
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
            loss$mz <- loss$mz - mass
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

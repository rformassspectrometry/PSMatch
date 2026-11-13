##' @title findBetterMatch
##'
##' Try to find a better peptide-spectrum match than the current PSM.
##'
##' This function is meant to try your luck at finding a better peptide-spectrum
##' match on the argument that [checkOverlap()] has returned false or that the
##' PSM's FDR is on the lower side. A peptide sequence is returned (with or
##' without modifications) for which the greatest number of matched peaks/a
##' better ndotproduct (see [MsCoreUtils::ndotproduct()] for more details) was
##' obtained.
##'
##' @param x a `Spectra()` object containing identifications.
##'
##' @param peptide a `character(1L)` pointing to the peptide sequence variable
##'   in `x`.
##'
##' @param modificationsList a `character()` with all modifications to try out.
##'   The modifications need to be present in the [modifications] dataset from
##'   the unimod database. Note that the more modifications are included, the
##'   more combinations are possible and the longer it takes to run.
##'
##' @param maxMods `integer(1L)` indicating the maximum modifications allowed at
##'   the same time on one sequence.
##'
##' @param improvedBy `character(1L)`: either `"maxScore"` or `"maxCount"`,
##'   defining how to choose the better match. `"maxScore"` keeps the match with
##'   the highest `ndotproduct` and `"maxCount"` keeps the one with the most
##'   matched peaks.
##'
##' @param convertToStyle
##'
##' @return A named `character()` vector of length equal to the number of
##'   spectra with identifications. The values are the best matches found using
##'   the `modificationsList` and the names are the indices for those sequences
##'   in the innitial Spectra object `x`.
##'
##' @author Guillaume Deflandre
##'
##' @examples
##' ## Load the data: a Spectra object with identifications
##' data(spBoekweg)
##' data(psmBoekweg)
##' head(psmBoekweg$pkey <- paste0(basename(psmBoekweg$filename),
##'                                sub("^.+scan=", "::", psmBoekweg$scannr)))
##' head(spBoekweg$pkey <- paste0(basename(spBoekweg$dataOrigin),
##'                               sub("^.+scan=", "::", spBoekweg$spectrumId)))
##' sp <- Spectra::joinSpectraData(spBoekweg, psmBoekweg, by.x = "pkey")
##'
##' ## Run with stock modifications (see PSMatch:::.ptmNames)
##' findBetterMatch(sp[30], "peptide", maxMods = 1)
##'
##' ## Run with custom modifications (need to be present in `modifications$Name`)
##' sp[30]$sequence <- getCanonicalSequence(sp[30]$peptide)
##' findBetterMatch(sp[30], "sequence", maxMods = 1)
##'
##' ## Get the best match based on the number of matched peaks
##' findBetterMatch(sp[30], "sequence", maxMods = 1, improvedBy = "maxCount")
findBetterMatch <- function(x,
                            peptide = "peptide",
                            modificationsList = NULL,
                            maxMods = 2,
                            improvedBy = c("maxScore", "maxCount"),
                            convertToStyle = "name") {
    ## Stop if x is not spectra object or if peptide sequence is not found:
    if (class(x) != "Spectra" | !("sequence" %in% Spectra::spectraVariables(x))) {
        stop(paste0("'x' must be a Spectra object and \nmust have an existing",
                    "'peptide' argument in Spectra::spectraVariables(x)"))
    }

    if (!requireNamespace("MsCoreUtils", quietly = TRUE))
            stop("Package 'MsCoreUtils' required for similarity calculation")

    improvedBy <- match.arg(improvedBy)

    if (is.null(modificationsList)) {
        modificationsList <- .ptmNames
    }

    ## Only keep identified spectra
    ids <- which(!is.na(x[[peptide]]))
    x <- x[ids, ]

    ## Get the canonical sequences on which to apply the function
    x$pepSequence <- getCanonicalSequence(x[[peptide]])

    ## Load data and select only a portion of it
    data(modifications)
    restrictedList <- modifications[modifications$Name %in% modificationsList, c("Name", "Site")]

    improvedSequences <- NULL
    for (i in seq_along(x)) {
        improvedSequences <- c(improvedSequences,
                               .findBetterMatch(x[i],
                                                restrictedList,
                                                maxMods,
                                                improvedBy = improvedBy))
    }

    ## Finally, convert the annotation style
    improvedSequences <- convertAnnotation(improvedSequences, convertToStyle = convertToStyle)
    names(improvedSequences) <- ids
    return(improvedSequences)
}

##' Helper function to `findBetterMatch()`
##'
##' @param x
##' @param pepSequences
##' @param restrictedList
##' @param maxMods
##' @return
##' @author GuillaumeDeflandre
##' @noRd
.findBetterMatch <- function(x,
                             restrictedList,
                             maxMods,
                             improvedBy = c("maxScore", "maxCounts")) {
    improvedBy <- match.arg(improvedBy)
    peptideSequences <- x$pepSequence

    for (ptm in unique(restrictedList$Name)) {
        sites <- restrictedList[restrictedList$Name == ptm, "Site"]
        ptm <- rep(ptm, length(sites))
        names(ptm) <- sites
        ptm <- ptm[!duplicated(names(ptm))]

        peptideSequences <- base::append(peptideSequences,
                                         addVariableModifications(peptideSequences,
                                                                  ptm,
                                                                  maxMods = maxMods))
        peptideSequences <- unique(peptideSequences)
    }

    ## Get observed spectrum
    observed <- Spectra::peaksData(x)

    ## Calculate theoretical sequence using labelFragments
    repeatedSpectra <- x[rep(seq_along(x), length(peptideSequences))]
    repeatedSpectra$sequence <- peptideSequences
    seqLabs <- suppressWarnings(
        labelFragments(repeatedSpectra, what = "mz")
    )

    similarities <- mapply(.calculateSimilarity,
                           seqLabs = seqLabs,
                           observed = observed)

    ## Select the improved sequence based on the choice If equal chances, take
    ## the first one (is more likely to be the most canonical sequence)
    ## Maybe keep all if names(improvedSequences) is set anyway ? TODO
    column <- if (improvedBy == "maxScore") 1L else 2L
    n <- which.max(similarities[column, ])[1]
    improvedSequence <- colnames(similarities)[n]
    return(improvedSequence)
}

##' Helper function: calculate similarity between observed and
##' theoretical spectrum using the ndotproduct
##'
##' @title
##' @param seqLabs
##' @param observed
##' @return
##' @author GuillaumeDeflandre
.calculateSimilarity <- function(seqLabs, observed) {
    obsMz <- observed[, "mz"]
    obsInt <- observed[, "intensity"]
    theoMz <- seqLabs

    ## Count matches
    nMatched <- sum(!is.na(theoMz))

    ## Prepare theoretical spectrum with matched intensities
    matchedIdx <- which(!is.na(theoMz))
    if (length(matchedIdx) == 0L) {
        return(list(similarity = 0, nMatched = 0L))
    }

    ## Create theoretical intensity vector (matched peaks get observed
    ## intensity)
    theoInt <- rep(0, length(obsMz))
    theoInt[matchedIdx] <- obsInt[matchedIdx]

    theoSpectrum <- matrix(c(theoMz, theoInt),
                        ncol = 2,
                        dimnames = list(c(), c("mz", "intensity")))

    ## Calculate cosine similarity using ndotproduct
    similarity <- MsCoreUtils::ndotproduct(observed, theoSpectrum)

    ## Store both similarities and number of matched peaks
    return(list(similarity = similarity, nMatched = nMatched))
}

## The stored restrictedList of modifications used in `findBetterMatch()`
.ptmNames <- c(
    # Phosphorylation and acidic modifications
    "Phospho",
    "Sulfo",
    "Carboxy",
    "Succinyl",
    "Malonyl",

    # Acetylation and acyl modifications
    "Acetyl",
    "Formyl",
    "Propionyl",
    "Butyryl",
    "Crotonyl",
    "Benzoyl",
    "Palmitoyl",
    "Myristoyl",
    "Octanoyl",
    "Lipoyl",
    "Biotin",

    # Methylation and related
    "Methyl",
    "Dimethyl",
    "Trimethyl",
    "Methylthio",

    # Lipid and isoprenoid modifications
    "Farnesyl",
    "GeranylGeranyl",

    # Glycosylation
    "Hex",
    "HexNAc",
    "Hex(1)HexNAc(1)",
    "NeuAc",
    "NeuGc",
    "ADP-Ribosyl",

    # Ubiquitin-like and remnant tags
    "LRGG",

    # Oxidation and redox
    "Oxidation",
    "Dioxidation",
    "Trioxidation",
    "Nitro",
    "Nitrosyl",
    "Dehydrated",
    "Sulfide",
    "Glutathione",
    "Cysteinyl",

    # Amide, deamidation and cyclisation
    "Amidated",
    "Deamidated",
    "Deamidated:18O(1)",
    "Carbamidomethyl",
    "Carbamyl",
    "Glu->pyro-Glu",
    "Gln->pyro-Glu",

    # Other common modifications
    "Phosphopantetheine",
    "Ammonia-loss",
    "Iodo",
    "Bromo"#,

    ## # Labelling and reagents (can't apply mods on extremities here..)
    ## "iTRAQ4plex",
    ## "iTRAQ8plex",
    ## "TMT6plex",
    ## "TMTpro",
    ## "Label:13C(6)",
    ## "Label:13C(6)15N(2)"
)

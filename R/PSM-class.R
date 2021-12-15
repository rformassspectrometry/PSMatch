##' @title A class for peptide-spectrum matches
##'
##' @aliases PSM-class PSM reducePSMs readPSMs psmVariables
##'
##' @name PSM
##'
##' @description
##'
##' The `PSM` class is a simple class to store and manipulate
##' peptide-spectrum matches. The class encapsulates PSM data as a
##' `DataFrame` (or more specifically a `DFrame`) with additional
##' lightweight metadata annotation.
##'
##' There are two types of `PSM` objects:
##'
##' - Objects with duplicated spectrum identifiers. This holds for
##'   multiple matches to the same spectrum, be it different peptide
##'   sequences or the same sequence with or without a
##'   post-translational modification. Such objects are typically
##'   created with the `PSM()` constructor starting from `mzIdentML`
##'   files.
##'
##' - Reduced objects where the spectrum identifiers (or any
##'   equivalent column) are unique keys within the PSM table. Matches
##'   to the same scan/spectrum are merged into a single PSM data
##'   row. Reduced `PSM` object are created with the `reducePSMs()`
##'   function. See examples below.
##'
##' Objects can be checked for their reduced state with the
##' `reduced()` function which returns `TRUE` for reduced instances,
##' `FALSE` when the spectrum identifiers are duplicated, or NA when
##' unknown. The flag can also be set explicitly with the
##' `reduced()<-` setter.
##'
##' @section Creating and using PSM objects:
##'
##' - The [PSM()] constructor uses parsers provided by the `mzR` or
##'   `mzID` packages to read the `mzIdentML` data. The vignette
##'   describes some apparent differences in their outputs. The
##'   constructor input is a character of one more multiple file
##'   names.
##'
##' - `PSM` objects can also be created from a `data.frame` object (or
##'   any variable that can be coerced into a [DataFrame].
##'
##' - Finally, [PSM()] can also take a `PSM` object as input, which
##'   leaves the PSM data as is and is used to set/update the PSM
##'   variables.
##'
##' - The constructor can also initialise variables (called *PSM
##'   variables*) needed for downstream processing, notably filtering
##'   (see [filterPSMs()]) and to generate a peptide-by-protein
##'   adjacency matrix (see [makeAdjacencyMatrix()]). These variables
##'   can be extracted with the [psmVariables()] function. They
##'   represent the columns in the PSM table that identify spectra,
##'   peptides, proteins, decoy peptides hit ranks and, optionally, a
##'   PSM score. The value of these variables will depend on the
##'   backend used to create the object, or left blank (i.e. encoded
##'   as `NA`) when building an object by hand from a `data.frame`. In
##'   such situation, they need to be passed explicitly by the user as
##'   arguments to [PSM()].
##'
##' - The `adjacencyMatrix()` accessor can be used to retrieve the
##'   binary sparse peptide-by-protein adjacency matrix from the PSM
##'   object. It also relies on PSM variables which thus need to be
##'   set beforehand. For more flexibility in the generation of the
##'   adjacency matrix (for non-binary matrices), use
##'   [makeAdjacencyMatrix()].
##'
##' @examples
##'
##' ## ---------------------------------
##' ## Example with a single mzid file
##' ## ---------------------------------
##'
##' f <- msdata::ident(full.names = TRUE, pattern = "TMT")
##' basename(f)
##'
##' ## mzR parser (default)
##' psm <- PSM(f)
##' psm
##'
##' ## PSM variables
##' psmVariables(psm)
##'
##' ## mzID parser
##' psm_mzid <- PSM(f, parser = "mzID")
##' psm_mzid
##'
##' ## different PSM variables
##' psmVariables(psm_mzid)
##'
##' ## Reducing the PSM data
##' (i <- which(duplicated(psm$spectrumID))[1:2])
##' (i <- which(psm$spectrumID %in% psm$spectrumID[i]))
##' psm2 <- psm[i, ]
##' reduced(psm2)
##'
##' ## Peptide sequence CIDRARHVEVQIFGDGKGRVVALGERDCSLQRR with
##' ## Carbamidomethyl modifications at positions 1 and 28.
##' DataFrame(psm2[, c("sequence", "spectrumID", "modName", "modLocation")])
##' reduced(psm2) <- FALSE
##' reduced(psm2)
##'
##' ## uses by default the spectrum PSM variable, as defined during
##' ## the construction - see psmVariables()
##' rpsm2 <- reducePSMs(psm2)
##' rpsm2
##' DataFrame(rpsm2[, c("sequence", "spectrumID", "modName", "modLocation")])
##' reduced(rpsm2)
##'
##' ## ---------------------------------
##' ## Multiple mzid files
##' ## ---------------------------------
##'
##' library(rpx)
##' PXD022816 <- PXDataset("PXD022816")
##' PXD022816
##'
##' (mzids <- pxget(PXD022816, grep("mzID", pxfiles(PXD022816))[1:2]))
##' psm <- PSM(mzids)
##' psm
##' psmVariables(psm)
##'
##' ## Here, spectrum identifiers are repeated accross files
##' psm[grep("scan=20000", psm$spectrumID), "spectrumFile"]
##'
##' ## Let's create a new primary identifier composed of the scan
##' ## number and the file name
##' psm$pkey <- paste(sub("^.+Task\\\\", "", psm$spectrumFile),
##'                   sub("^.+scan=", "", psm$spectrumID),
##'                   sep = "::")
##' head(psm$pkey)
##'
##' ## the PSM is not reduced
##' reduced(psm, "pkey")
##' DataFrame(psm[6:7, ])
##'
##' ## same sequence, same spectrumID, same file
##' psm$sequence[6:7]
##' psm$pkey[6:7]
##'
##' ## different modification locations
##' psm$modLocation[6:7]
##'
##' ## here, we need to *explicitly* set pkey to reduce
##' rpsm <- reducePSMs(psm, psm$pkey)
##' rpsm
##' reduced(rpsm, "pkey")
##'
##' ## the two rows are now merged into a single one; the distinct
##' ## modification locations are preserved.
##' (i <- which(rpsm$pkey == "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::12894"))
##' DataFrame(rpsm[i, c("sequence", "pkey", "modName", "modLocation")])
##'
##' ## ---------------------------------
##' ## PSM from a data.frame
##' ## ---------------------------------
##'
##' psmdf <- data.frame(spectrum = paste0("sp", 1:10),
##'                     sequence = replicate(10,
##'                                          paste(sample(getAminoAcids()[-1, "AA"], 10),
##'                                                collapse = "")),
##'                     protein = sample(paste0("Prot", LETTERS[1:7]), 10,
##'                                      replace = TRUE),
##'                     decoy = rep(FALSE, 10),
##'                     rank = rep(1, 10),
##'                     score = runif(10))
##' psmdf
##'
##' psm <- PSM(psmdf)
##' psm
##' psmVariables(psm)
##'
##' ## no PSM variables set
##' try(adjacencyMatrix(psm))
##'
##' ## set PSM variables
##' psm <- PSM(psm, spectrum = "spectrum", peptide = "sequence",
##'            protein = "protein", decoy = "decoy", rank = "rank")
##' psm
##' psmVariables(psm)
##'
##' adjacencyMatrix(psm)
NULL


##' @name PSM-class
##'
##' @docType class
##'
##' @exportClass PSM
##'
##' @importClassesFrom S4Vectors DataFrame DFrame
##'
##' @noRd
setClass("PSM",
         contains = "DFrame")

showDetails <- function(object) {
    .psmVariables <- psmVariables(object)
    n_sp <- length(unique(object[[.psmVariables["spectrum"]]]))
    n_pe <- length(unique(object[[.psmVariables["peptide"]]]))
    n_pr <- length(unique(object[[.psmVariables["protein"]]]))
    tabDecoy <- object[[.psmVariables["decoy"]]]
    cat("Spectra:", n_sp, "unique\n")
    cat("  db: ", sum(!tabDecoy), " target, ",
        sum(tabDecoy), " decoy\n", sep = "")
    tabRank <- table(object[[.psmVariables["rank"]]])
    cat("  ranks:", paste0(names(tabRank),":", tabRank), "\n")
    cat("Peptides: ")
    mlt <- tapply(object[[psmVariables(object)["protein"]]],
                  object[[psmVariables(object)["peptide"]]],
                  function(xx) length(unique(xx)) > 1)
    cat(sum(!mlt), "unique,", sum(mlt), "multiple\n")
    cat("Proteins:", n_pr, "\n")
    invisible(NULL)
}

setMethod("show", "PSM",
          function(object) {
              cl <- classNameForDisplay(object)
              if (isTRUE(reduced(object)))
                  cl <- paste("Reduced", cl)
              cat(cl, "with", nrow(object), "rows and",
                  ncol(object), "columns.\n")
              if (!any(is.na(psmVariables(object))) & !isTRUE(reduced(object)) & nrow(object) > 0)
                  showDetails(object)
              if (ncol(object) <= 4)
                  cat("names(", ncol(object), "): ",
                      paste(names(object), collapse = " "), "\n",
                      sep = "")
              else
                  cat("names(", ncol(object), "): ",
                      paste(names(object)[1:2], collapse = " "),
                      " ... ", paste(names(object)[(ncol(object)-1):ncol(object)],
                                   collapse = " "), "\n",
                      sep = "")
              invisible(NULL)
          })




##' @param x `character()` of mzid file names, an instance of class
##'     `PSM`, or a `data.frame`.
##'
##' @param spectrum `character(1)` variable name that defines a
##'     spectrum in the PSM data. Default are `"spectrumID"` (mzR
##'     parser) or `"spectrumid"` (mzID parser). It is also used to
##'     calculate the reduced state.
##'
##' @param peptide `character(1)` variable name that defines a peptide
##'     in the PSM data. Detaults are `"spequence"` (mzR parser) or
##'     `"pepSeq"` (mzID parser).
##'
##' @param protein `character(1)` variable name that defines a protein
##'     in the PSM data. Detaults are `"DatabaseAccess"` (mzR parser)
##'     or `"accession"` (mzID parser).
##'
##' @param decoy `character(1)` variable name that defines a decoy hit
##'     in the PSM data. Detaults are `"isDecoy"` (mzR parser) or
##'     `"isdecoy"` (mzID parser).
##'
##' @param rank `character(1)` variable name that defines the rank of
##'     the peptide spectrum match in the PSM data. Default is `"rank"`.
##'
##' @param score `character(1)` variable name that defines the PSM
##'     score. This value isn't set by default as it depends on the
##'     search engine and application. Default is `NA`.
##'
##' @param parser `character(1)` defining the parser to be used to
##'     read the `mzIdentML` files. One of `"mzR"` (default) or
##'     `"mzID"`.
##'
##' @param BPPARAM an object inheriting from [BiocParallelParam] to
##'     control parallel processing. The default value is
##'     `SerialParam()` to read files in series.
##'
##' @return `PSM()` returns a `PSM` object.
##'
##' @export PSM
##'
##' @aliases PSM PSM,character PSM,data.frame PSM,PSM
##'
##' @import S4Vectors BiocParallel
##'
##' @name PSM
PSM <- function(x,
                spectrum = NA,
                peptide = NA,
                protein = NA,
                decoy = NA,
                rank = NA,
                score = NA,
                parser = c("mzR", "mzID"),
                BPPARAM = SerialParam()) {
    if (is.character(x)) {
        if (!all(flex <- file.exists(x)))
            stop(paste(x[!flex], collapse = ", "), " not found.")
        parser <- match.arg(parser)
        if (parser == "mzR") {
            psm <- readPSMsMzR(x, BPPARAM)
            ## Default PSM variables for mzR parser
            .psmVariables <- c(spectrum = "spectrumID",
                               peptide = "sequence",
                               protein = "DatabaseAccess",
                               decoy = "isDecoy",
                               rank = "rank",
                               score = NA_character_)
        } else {
            psm <- readPSMsMzID(x, BPPARAM)
            ## Default PSM variables for mzID parser
            .psmVariables <- c(spectrum = "spectrumid",
                               peptide = "pepseq",
                               protein = "accession",
                               decoy = "isdecoy",
                               rank = "rank",
                               score = NA_character_)
        }
    } else if (is.data.frame(x)) {
        psm <- as(DataFrame(x), "PSM")
        .psmVariables <- c(spectrum = NA_character_,
                           peptide = NA_character_,
                           protein = NA_character_,
                           decoy = NA_character_,
                           rank = NA_character_,
                           score = NA_character_)
    } else {
        stopifnot(inherits(x, "PSM"))
        .psmVariables <- psmVariables(x)
        psm <- x ## use as is
    }
    ## Update PSM variables based on contstructor inputs
    if (spectrum %in% names(psm))
        .psmVariables["spectrum"] <- spectrum
    if (peptide %in% names(psm))
        .psmVariables["peptide"] <- peptide
    if (protein %in% names(psm))
        .psmVariables["protein"] <- protein
    if (decoy %in% names(psm))
        .psmVariables["decoy"] <- decoy
    if (rank %in% names(psm))
        .psmVariables["rank"] <- rank
    if (score %in% names(psm))
        .psmVariables["score"] <- score
    metadata(psm)$variables <- .psmVariables
    psm
}

##' @name PSM
##'
##' @export
reduced <- function(object, spectrum = psmVariables(object)["spectrum"]) {
    if (!missing(spectrum))
        return(!is.list(object[[spectrum]]) & !anyDuplicated(object[[spectrum]]))
    if (is.null(metadata(object)[["reduced"]]))
        return(NA)
    else metadata(object)[["reduced"]]
}

##' @param value new value to be passed to setter.
##'
##' @name PSM
##'
##' @export
"reduced<-" <- function(object, value) {
    stopifnot(is.logical(value))
    metadata(object)[["reduced"]] <- value
    object
}

##' @param object An instance of class `PSM`.
##'
##' @param which `character()` with the PSM variable name to
##'     retrieve. If `"all"` (default), all named variables are
##'     returned. See [PSM()] for valid PSM variables.
##'
##' @name PSM
##'
##' @export
psmVariables <- function(object, which = "all") {
    .psmVariables <- metadata(object)[["variables"]]
    if (length(which) == 1 && which == "all")
        return(.psmVariables)
    stopifnot(which %in% names(.psmVariables))
    .psmVariables[which]
}

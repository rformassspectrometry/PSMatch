##' @title A class for peptide-spectrum matches
##'
##' @aliases PSM-class PSM reducePSMs readPSMs
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
##' - Reduced object where the spectrum identifiers (or any equivalent
##'   column - see example below) are unique keys within the PSM
##'   table. Matches to the same scan/spectrum in each raw data file
##'   are merged into a single data row. Reduced `PSM` object are
##'   created with the `reducePSMs()` function. See examples below.
##'
##' Objects can be checked for their reduced state with the
##' `reduced()` function which returns TRUE for reduced instances,
##' FALSE when the spectrum identifiers are duplicated, or NA when
##' unknown. The flag can also be set explicitly with the
##' `reduced()<-` setter.
##'
##' The constructor uses parsers provided by the `mzR` or `mzID`
##' packages to read the `mzIdentML` data. See the vignette for some
##' apparent differences in their outputs.
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
##' ## mzID parser
##' PSM(f, parser = "mzID")
##'
##' i <- which(duplicated(psm$spectrumID))[1:2]
##' i <- which(psm$spectrumID %in% psm$spectrumID[i])
##' psm2 <- psm[i, ]
##' reduced(psm2)
##'
##' ## Peptide sequence CIDRARHVEVQIFGDGKGRVVALGERDCSLQRR with
##' ## Carbamidomethyl modifications at positions 1 and 28.
##' DataFrame(psm2[, c("sequence", "spectrumID", "modName", "modLocation")])
##' reduced(psm2) <- FALSE
##'
##' rpsm2 <- reducePSMs(psm2, psm2$spectrumID)
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
##' rpsm <- reducePSMs(psm, psm$pkey)
##' rpsm
##' reduced(rpsm, "pkey")
##'
##' ## the two rows are now merged into a single one; the distinct
##' ## modification locations are preserved.
##' (i <- which(rpsm$pkey == "QEP2LC6_HeLa_50ng_251120_01-calib.mzML::12894"))
##' DataFrame(rpsm[i, c("sequence", "pkey", "modName", "modLocation")])
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

setMethod("show", "PSM",
          function(object) {
              cl <- classNameForDisplay(object)
              if (isTRUE(reduced(object)))
                  cl <- paste("Reduced", cl)
              cat(cl, "with", nrow(object), "rows and",
                  ncol(object), "columns.\n")
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


##' @param x `character()` of `mzid` file names or an instance of
##'     class `PSM`.
##'
##' @param spectrum `character(1)` that defines a spectrum in the PSM
##'     data. Typically `"spectrumID"`.
##'
##' @param peptide `character(1)` that defines a peptide in the PSM
##'     data. Typically `"spequence"`.
##'
##' @param protein `character(1)` that defines a protein in the PSM
##'     data. Typically `"DatabaseAccess"`.
##'
##' @param decoy `character(1)` that defines a decoy hit in the PSM
##'     data. Default is `"isDecoy"`. This variable is named
##'     `"isdecoy"` if you use the `mzID` parser.
##'
##' @param rank `character(1)` that defines the rank of the petide
##'     spetrum match in the PSM data. Default is `"rank"`.
##'
##' @param parser `character(1)` defining the parser to be used to
##'     read the `mzIdentML` files. One of `"mzR"` (default) or
##'     `"mzID"`.
##'
##' @param BPPARAM an object inheriting from [BiocParallelParam] to
##'     control parallel processing. The default value is
##'     `SerialParam()` to read files in series.
##'
##' @return `PSM()` returns a `PSM` object containing the
##'     peptide-spectrum matches stored in `files`.
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
                decoy = "isDecoy",
                rank = "rank",
                parser = c("mzR", "mzID"),
                BPPARAM = SerialParam()) {
    if (is.character(x)) {
        if (!all(flex <- file.exists(x)))
            stop(paste(x[!flex], collapse = ", "), " not found.")
        parser <- match.arg(parser)
        if (parser == "mzR") psm <- readPSMsMzR(x, BPPARAM)
        else psm <- readPSMsMzID(x, BPPARAM)
    } else if (is.data.frame(x)) {
        ## TODO - constructor for data.frame
    } else {
        stopifnot(inherits(x, "PSM"))
        psm <- x
    }
    ## Store usefull variables if provided as part of the
    ## constructor. There could also be an interactive function that
    ## asks the user to set these.
    psmVariables <- c(spectrum = NA_character_,
                      peptide = NA_character_,
                      protein = NA_character_,
                      decoy = NA_character_,
                      rank = NA_character_)
    metadata(psm)$variables <- psmVariables
    if (spectrum %in% names(psm))
        metadata(psm)$variables["spectrum"] <- spectrum
    if (peptide %in% names(psm))
        metadata(psm)$variables["peptide"] <- peptide
    if (protein %in% names(psm))
        metadata(psm)$variables["protein"] <- protein
    if (decoy %in% names(psm))
        metadata(psm)$variables["decoy"] <- decoy
    if (rank %in% names(psm))
        metadata(psm)$variables["rank"] <- rank
    psm
}



##' @param spectrumID `character(1)` with the column name referring to
##'     the spectrum identifier or any other unique key for the `PSM`
##'     object `x`. This parameter can be set to calculate the reduced
##'     state explicitly by by-passing the `reduced` metadata flag.
##'
##' @param x An instance of class `PSM`.
##'
##' @name PSM
##'
##' @export
reduced <- function(x, spectrumID) {
    if (!missing(spectrumID))
        return(!is.list(x[[spectrumID]]) & !anyDuplicated(x[[spectrumID]]))
    if (is.null(metadata(x)[["reduced"]]))
        return(NA)
    else metadata(x)[["reduced"]]
}

##' @param value `logical(1)` to set the reduced state of the `PSM`
##'     object.
##'
##' @name PSM
##'
##' @export
"reduced<-" <- function(x, value) {
    stopifnot(is.logical(value))
    metadata(x)[["reduced"]] <- value
    x
}

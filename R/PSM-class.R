##' @title A class for peptide-spectrum matches results
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


##' @param files `character()` of `mzid` file names.
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
##' @import S4Vectors BiocParallel
##'
##' @name PSM
PSM <- function(files,
                parser = c("mzR", "mzID"),
                BPPARAM = SerialParam()) {
    if (!all(flex <- file.exists(files)))
        stop(paste(files[!flex], collapse = ", "), " not found.")
    parser <- match.arg(parser)
    if (parser == "mzR") readPSMsMzR(files, BPPARAM)
    else readPSMsMzID(files, BPPARAM)
}

##' @param spectrumID `character(1)` with the column name referring to
##'     the spectrum identifier or any other unique key for the `PSM`
##'     object `x`.
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

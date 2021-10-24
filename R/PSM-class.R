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
##' `reduced()` function which returns `TRUE` for reduced instances,
##' `FALSE` when the spectrum identifiers are duplicated, or NA when
##' unknown. The flag can also be set explicitly with the
##' `reduced()<-` setter.
##'
##' - The constructor can also be initialise some variables needed for
##'   downstream processing, notably filtering (See
##'   [filterPSMs()]). These variables can be extracted with the
##'   [psmVariables()] function. They represent the columns in the PSM
##'   table that identify spectra, peptides, proteins, decoy peptides
##'   and hit ranks.
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
##' ## PSM variables
##' psmVariables(psm)
##'
##' ## mzID parser
##' psm_mzid <- PSM(f, parser = "mzID")
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
    cat("Spectra: ")
    cat(sum(!tabDecoy), " target, ",
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




##' @param x `character()` of `mzid` file names or an instance of
##'     class `PSM`.
##'
##' @param spectrum `character(1)` that defines a spectrum in the PSM
##'     data. Default are `"spectrumID"` (mzR parser) or
##'     `"spectrumid"` (mzID parser). It is also used to calculate the
##'     reduced state.
##'
##' @param peptide `character(1)` that defines a peptide in the PSM
##'     data. Detaults are `"spequence"` (mzR parser) or `"pepSeq"`
##'     (mzID parser).
##'
##' @param protein `character(1)` that defines a protein in the PSM
##'     data. Detaults are `"DatabaseAccess"` (mzR parser) or
##'     `"accession"` (mzID parser).
##'
##' @param decoy `character(1)` that defines a decoy hit in the PSM
##'     data. Detaults are `"isDecoy"` (mzR parser) or `"isdecoy"`
##'     (mzID parser).
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
                decoy = NA,
                rank = NA,
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
                               rank = "rank")
        } else {
            psm <- readPSMsMzID(x, BPPARAM)
            ## Default PSM variables for mzID parser
            .psmVariables <- c(spectrum = "spectrumid",
                               peptide = "pepseq",
                               protein = "accession",
                               decoy = "isdecoy",
                               rank = "rank")
        }
    } else if (is.data.frame(x)) {
        ## TODO - constructor for data.frame
    } else {
        stopifnot(inherits(x, "PSM"))
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
        psmVariables["rank"] <- rank
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

##' @name PSM
##'
##' @export
"psmVariables<-" <- function(object, value) {
    value <- as.character(value)
    metadata(object)[["variables"]] <- value
    object
}

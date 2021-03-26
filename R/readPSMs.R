##' Reads a set of `mzId` files containing PSMs and returns the PSMs
##' as a `DataFrame`.
##'
##' This function uses parsers provided by the `mzR` or `mzID`
##' packages to read the `mzIdentML` data. See the vignette for some
##' apparent differences in their outputs.
##'
##' @title Import peptide-spectrum matches
##'
##' @param files A `character` of `mzid` files.
##'
##' @param parser `character(1)` defining the parser to be used to
##'     read the `mzIdentML` files. One of `"mzR"` (default) or
##'     `"mzID"`.
##'
##' @param BPPARAM an object inheriting from [BiocParallelParam] to
##'     control parallel processing. The default value is
##'     `SerialParam()` to read files in series.
##'
##' @return A `DataFrame` containing the PSMs stored in the `mzId`
##'     files.
##'
##' @author Laurent Gatto
##'
##' @export
##'
##' @import S4Vectors BiocParallel
##'
##' @examples
##' f <- msdata::ident(full.names = TRUE, pattern = "TMT")
##' basename(f)
##'
##' ## mzR parser (default)
##' readPSMs(f)
##'
##' ## mzID parser
##' readPSMs(f, parser = "mzID")
readPSMs <- function(files,
                     parser = c("mzR", "mzID"),
                     BPPARAM = SerialParam()) {
    if (!all(flex <- file.exists(files)))
        stop(paste(files[!flex], collapse = ", "), " not found.")
    parser <- match.arg(parser)
    if (parser == "mzR") readPSMsMzR(files, BPPARAM)
    else readPSMsMzID(files, BPPARAM)
}


##' @importFrom methods as
readPSMsMzR <- function(files, BPPARAM) {
    stopifnot(requireNamespace("mzR"))
    iddf <- bplapply(files,
                     function(f) as_data_frame(mzR::openIDfile(f)),
                     BPPARAM = BPPARAM)
    iddf <- do.call(rbind, iddf)
    iddf <- as(iddf, "DataFrame")
    as(iddf, "PSM")
}


## mzID automatically reads data in parallel and uses all the cores,
## without giving the user any control. This function doesn't make use
## of this by passing only one file at a time to the constructor.
readPSMsMzID <- function(files, BPPARAM) {
    stopifnot(requireNamespace("mzID"))
    iddf <- bplapply(files,
                     function(f) mzID::flatten(mzID::mzID(f)),
                     BPPARAM = BPPARAM)
    iddf <- do.call(rbind, iddf)
    iddf <- as(iddf, "DataFrame")
    as(iddf, "PSM")
}

##' @export PSM
PSM <- readPSMs

##' Reads as set of `mzId` files containing PSMs and returns the PSMs
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
##' @return A `DataFrame` containing the PSMs stored in the `mzId`
##'     files.
##' 
##' @author Laurent Gatto
##'
##' @export
##'
##' @import S4Vectors
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
readPSMs <- function(files, parser = c("mzR", "mzID")) {
    if (!all(flex <- file.exists(files)))
        stop(paste(files[!flex], collapse = ", "), " not found.")    
    parser <- match.arg(parser)
    if (parser == "mzR") readPSMsMzR(files)
    else readPSMsMzID(files)
}


##' @importFrom methods as
readPSMsMzR <- function(files) {
    stopifnot(requireNamespace("mzR"))
    if (length(files) == 1) {
        id <- mzR::openIDfile(files)
        iddf <- as(id, "data.frame")
    } else {
        iddf <- lapply(files,
                       function(f) as(mzR::openIDfile(f), "data.frame"))
        iddf <- do.call(rbind, iddf)
    }
    as(iddf, "DataFrame")
}


readPSMsMzID <- function(files) {
    stopifnot(requireNamespace("mzID"))
    if (length(files) == 1) {
        iddf <- mzID::flatten(mzID::mzID(files))
    } else {
        iddf <- lapply(files,
                       function(f) mzID::flatten(mzID::mzID(files)))
        iddf <- do.call(rbind, iddf)
    }    
    as(iddf, "DataFrame")
}

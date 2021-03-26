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
## without giving the user any control. This function avoids this by
## passing only one file at a time to the constructor.
readPSMsMzID <- function(files, BPPARAM) {
    stopifnot(requireNamespace("mzID"))
    iddf <- bplapply(files,
                     function(f) mzID::flatten(mzID::mzID(f)),
                     BPPARAM = BPPARAM)
    iddf <- do.call(rbind, iddf)
    iddf <- as(iddf, "DataFrame")
    as(iddf, "PSM")
}


## for backward compatibility
##' @export
readPSMs <- PSM

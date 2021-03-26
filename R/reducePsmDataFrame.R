## add a simplify that will take unique when multiple identical
## elements

##' This function will reduce a DataFrame of peptide-spectrum matches,
##' as stored as a PSM object. Scans, as uniquely defined by their
##' spectrum identifiers, can match multiple peptides, be it different
##' peptide sequences or the same sequence with or without a
##' post-translational modification.
##'
##' @title Reduce a PSM DataFrame
##'
##' @param x
##'
##' @param k
##'
##' @return
##'
##' @author Laurent Gatto
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
##' i <- anyDuplicated(psm$spectrumID)
##' i <- which(psm$spectrumID %in% psm$spectrumID[i])
##' psm2 <- psm[i, ]
##'
##' ## Peptide sequence CIDRARHVEVQIFGDGKGRVVALGERDCSLQRR with
##' ## Carbamidomethyl modifications at positions 1 and 28.
##' as.data.frame(psm2)
reducePsmDataFrame <- function(x, k) {
    if (missing(k))
        stop("Argument k is missing")
    x <- QFeatures::reduceDataFrame(x, k)
    n <- ncol(x)
    for (i in seq_along(x)) {
        .x <- x[[i]]
        class_x <- class(.x)
        .x <- sapply(.x, unique)
        if (is.list(.x))
            .x <- as(.x, class_x)
        else .x <- unname(.x)
        x[[i]] <- .x
    }
    metadata(x)[["reduced"]] <- TRUE
    as(x, "PSM")
}

##' A function to filter out PSMs matching to the decoy database, of
##' rank greater than one and matching non-proteotypic peptides.
##'
##' The PSMs should be stored in a `DataFrame` (or `data.frame` or
##' `tibble`) such as those produced by [readPSMs()]. Note that
##' this function should be used before reducing the PSM table.
##'
##' @title Filter out unreliable PSMs.
##' 
##' @param x A `DataFrame` containing PSMs.
##' 
##' @param decoy The column name defining whether entries match the
##'     decoy database. Default is `"isDecoy"`. The column should be a
##'     `logical` and only PSMs holding a `FALSE` are
##'     retained. Ignored is set to `NULL`.
##' 
##' @param rank The column name holding the rank of the PSM. Default
##'     is `"rank"`. This column should be a `numeric` and only PSMs
##'     having rank equal to 1 are retained. Ignored is set to `NULL`.
##' 
##' @param accession The column name holding the protein (groups)
##'     accession. Default is `"DatabaseAccess"`. Ignored is set to
##'     `NULL`.
##' 
##' @param spectrumID The name of the spectrum identifier
##'     column. Default is `spectrumID`.
##' 
##' @param verbose `logical(1)` setting the verbosity flag.
##' 
##' @return A new `DataFrame` with filtered out peptides and with the
##'     same columns as the input `x`.
##' 
##' @author Laurent Gatto
##'
##' @export
##'
##' @examples
##' f <- msdata::ident(full.names = TRUE, pattern = "TMT")
##' basename(f)
##' id <- readPSMs(f)
##' filterPSMs(id)
filterPSMs <- function(x,
                       decoy = "isDecoy",
                       rank = "rank",
                       accession = "DatabaseAccess",
                       spectrumID = "spectrumID",
                       verbose = TRUE) {
    n0 <- nrow(x)
    if (verbose)
        message("Starting with ", n0, " PSMs:")
    if (!is.null(decoy)) {
        x <- x[!x[, decoy], ]
        n1 <- nrow(x)
        if (verbose)
            message(" removed ", n0 - n1, " decoy hits")
        n0 <- n1
    }
    if (!is.null(rank)) {
        x <- x[x[, rank] == 1, ]
        n2 <- nrow(x)
        if (verbose)
            message(" removed ", n0 - n2, " PSMs with rank > 1")
        n0 <- n2
    }
    if (!is.null(accession)) {
        mlt <- tapply(x[, accession],
                      x[, spectrumID],
                      function(xx) length(unique(xx)) > 1)
        mlt <- names(which(mlt))        
        x <- x[!x$spectrumID %in% mlt, ]
        n3 <- nrow(x)
        if (verbose)
            message(" removed ", n0 - n3, " non-proteotypic peptides")
    }
    if (verbose)
        message(nrow(x), " PSMs left.")
    x
}

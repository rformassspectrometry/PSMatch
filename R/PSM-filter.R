##' @title Filter out unreliable PSMs.
##'
##' @description
##'
##' Functions to filter out PSMs matching. The PSMs should be stored
##' in a `PSM` such as those produced by [PSM()].
##'
##' @param x An instance of class `PSM`.
##'
##' @param decoy `character(1)` with the column name specifying
##'     whether entries match the decoy database or not. Default is
##'     the `decoy` PSM variable as defined in [psmVariables()]. The
##'     column should be a `logical` and only PSMs holding a `FALSE`
##'     are retained. Filtering is ignored if set to `NULL` or `NA`.
##'
##' @param rank `character(1)` with the column name holding the rank
##'     of the PSM. Default is the `rank` PSM variable as defined in
##'     [psmVariables()]. This column should be a `numeric` and only
##'     PSMs having rank equal to 1 are retained. Filtering is ignored
##'     if set to `NULL` or `NA`.
##'
##' @param protein `character(1)` with the column name holding the
##'     protein (groups) protein. Default is the `protein` PSM
##'     variable as defined in [psmVariables()]. Filtering is ignored
##'     if set to `NULL` or `NA`.
##'
##' @param spectrum `character(1)` with the name of the spectrum
##'     identifier column. Default is the `spectrum` PSM variable as
##'     defined in [psmVariables()]. Filtering is ignored if set to
##'     `NULL` or `NA`.
##'
##' @param peptide `character(1)` with the name of the peptide
##'     identifier column. Default is the `peptide` PSM variable as
##'     defined in [psmVariables()]. Filtering is ignored if set to
##'     `NULL` or `NA`.
##'
##' @param verbose `logical(1)` setting the verbosity flag.
##'
##' @return A new filtered `PSM` object with the same columns as the
##'     input `x`.
##'
##' @author Laurent Gatto
##'
##' @name filterPSMs
##'
##' @export
##'
##' @examples
##' f <- msdata::ident(full.names = TRUE, pattern = "TMT")
##' basename(f)
##' id <- PSM(f)
##' filterPSMs(id)
filterPSMs <- function(x,
                       decoy = psmVariables(x)["decoy"],
                       rank = psmVariables(x)["rank"],
                       protein = psmVariables(x)["protein"],
                       spectrum = psmVariables(x)["spectrum"],
                       peptide = psmVariables(x)["peptide"],
                       verbose = TRUE) {
    n0 <- nrow(x)
    if (verbose)
        message("Starting with ", n0, " PSMs:")
    x <- filterPsmDecoy(x, decoy = decoy,
                        verbose = verbose)

    x <- filterPsmRank(x, rank = rank,
                       verbose = verbose)

    x <- filterPsmShared(x,
                         protein = protein,
                         peptide = peptide,
                         verbose = verbose)
    if (verbose)
        message(nrow(x), " PSMs left.")
    x
}


##' @description
##'
##' - `filterPsmDecoy()` filters out decoy PSMs, i.e. those annotated
##'    as `isDecoy`.
##'
##' @name filterPSMs
##'
##' @export
filterPsmDecoy <- function(x,
                           decoy = psmVariables(x)["decoy"],
                           verbose = TRUE) {
    if (is.null(decoy) || is.na(decoy))
        return(x)
    n0 <- nrow(x)
    x <- x[!x[, decoy], ]
    n1 <- nrow(x)
    if (verbose)
        message("Removed ", n0 - n1, " decoy hits.")
    x
}

##' @description
##'
##' - `filterPsmRank()` filters out PSMs of rank > 1.
##'
##' @name filterPSMs
##'
##' @export
filterPsmRank <- function(x,
                          rank = psmVariables(x)["rank"],
                          verbose = TRUE) {
    if (is.null(rank) || is.na(rank))
        return(x)
    n0 <- nrow(x)
    x <- x[x[, rank] == 1, ]
    n1 <- nrow(x)
    if (verbose)
        message("Removed ", n0 - n1, " PSMs with rank > 1.")
    x
}


##' @description
##'
##' - `filterPsmShared()` filters out shared PSMs, i.e. those that
##'    match multiple proteins.
##'
##' @name filterPSMs
##'
##' @export
filterPsmShared <- function(x,
                            protein = psmVariables(x)["protein"],
                            peptide = psmVariables(x)["peptide"],
                            verbose = TRUE) {
    if (is.null(protein) || is.null(peptide) || (is.na(protein) | is.na(peptide)))
        return(x)
    n0 <- nrow(x)
    mlt <- tapply(x[, protein],
                  x[, peptide],
                  function(xx) length(unique(xx)) > 1)
    mlt <- names(which(mlt))
    x <- x[!x[[peptide]] %in% mlt, ]
    n1 <- nrow(x)
    if (verbose)
        message("Removed ", n0 - n1, " shared peptides.")
    x
}

##' @description
##'
##' - `filterPsmFdr()` filters out PSMs based on their FDR.
##'
##' @name filterPSMs
##'
##' @param fdr `character(1)` variable name that defines that defines
##'     the spectrum FDR (or any similar/relevant metric that can be
##'     used for filtering). This value isn't set by default as it
##'     depends on the search engine and application. Default is `NA`.
##'
##' @param FDR `numeric(1)` to be used to filter based on the `fdr`
##'     variable. Default is 0.05.
##'
##' @export
filterPsmFdr <- function(x,
                         FDR = 0.05,
                         fdr = psmVariables(x)["fdr"],
                         verbose = TRUE) {
    if (is.null(fdr) || is.na(fdr))
        return(x)
    n0 <- nrow(x)
    x <- x[x[, fdr] < FDR, ]
    n1 <- nrow(x)
    if (verbose)
        message("Removed ", n0 - n1, " PSMs with FDR > ", FDR, ".") 
    x
}

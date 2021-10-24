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
##'     are retained. Filtering is ignored if set to `NULL`.
##'
##' @param rank `character(1)` with the column name holding the rank
##'     of the PSM. Default is the `rank` PSM variable as defined in
##'     [psmVariables()]. This column should be a `numeric` and only
##'     PSMs having rank equal to 1 are retained. Filtering is ignored
##'     if set to `NULL`.
##'
##' @param protein `character(1)` with the column name holding the
##'     protein (groups) protein. Default is the `protein` PSM
##'     variable as defined in [psmVariables()]. Filtering is ignored
##'     if set to `NULL`.
##'
##' @param spectrum `character(1)` with the name of the spectrum
##'     identifier column. Default is the `spectrum` PSM variable as
##'     defined in [psmVariables()]. Filtering is ignored if set to
##'     `NULL`.
##'
##' @param peptide `character(1)` with the name of the peptide
##'     identifier column. Default is the `peptide` PSM variable as
##'     defined in [psmVariables()]. Filtering is ignored if set to
##'     `NULL`.
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
    if (!is.null(decoy))
        x <- filterPsmDecoy(x, decoy = decoy,
                            verbose = verbose)

    if (!is.null(rank))
        x <- filterPsmRank(x, rank = rank,
                           verbose = verbose)

    if (!is.null(protein))
        x <- filterPsmNonProteotypic(x,
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
    if (is.null(decoy))
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
    if (is.null(rank))
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
##' - `filterPsmNonProteotypic()` filters out non-proteotypic PSMs,
##'    i.e. those that match multiple proteins.
##'
##' @name filterPSMs
##'
##' @export
filterPsmNonProteotypic <- function(x,
                                    protein = psmVariables(x)["protein"],
                                    peptide = psmVariables(x)["peptide"],
                                    verbose = TRUE) {
    if (is.null(protein) | is.null(peptide))
        return(x)
    n0 <- nrow(x)
    mlt <- tapply(x[, protein],
                  x[, peptide],
                  function(xx) length(unique(xx)) > 1)
    mlt <- names(which(mlt))
    x <- x[!x[[peptide]] %in% mlt, ]
    n1 <- nrow(x)
    if (verbose)
        message("Removed ", n0 - n1, " non-proteotypic peptides.")
    x
}

## filterPsmUniqueSeq <- function(x,
##                                sequence = "sequence",
##                                verbose = TRUE) {}


##' @description
##'
##' - `filterPsmBestScore()` filters out PSMs but those with he
##'    highest score when matching the same spectrum identifier.
##'
##' @param score `character(1)` with the name of a score column to be
##'     used. Filtering is ignored if set to `NULL`.
##'
##' @name filterPSMs
filterPsmBestScore <- function(x,
                               score,
                               spectrum = psmVariables(x)["spectrum"],
                               verbose = TRUE) {
    if (missing(score))
        stop("Please provide a score variable.")
    if (!score %in% names(x))
        stop("Please provide a valid score.")
    if (is.null(score) | is.null(spectrum))
        return(x)
    stop("TODO")
}



##' @description
##'
##' - `filterPsmMods()` filters out PSMs that contain any PTM, as
##'    defined by non-missing (i.e. non-`NA`) modification name.
##'
##' @param mod `character(1)` with the column name holding the
##'     modification name. Default is `"modName"`. Filtering is
##'     ignored if set to `NULL`.
##'
##' @name filterPSMs
##'
##' @export
filterPsmMods <- function(x,
                          mod = "modName",
                          verbose = TRUE) {
    if (is.null(mod))
        return(x)
    n0 <- nrow(x)
    x <- x[is.na(x[, mod]), ]
    n1 <- nrow(x)
    if (verbose)
        message("Removed ", n0 - n1, " modified peptides.")
    x
}

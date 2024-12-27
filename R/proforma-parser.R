#' ProForma parser
#'
#' This helper functions provide a basic implementation of the PSI ProForma
#' notation.
#'
#' @author Sebastian Gibb <mail@sebastiangibb.de>
#'
#' @references http://www.psidev.info/proforma
#' @name proforma-parser

#' Remove all modifications
#'
#' @name proforma-parser
#' @param x `character`, ProForma sequence.
#' @return `character`, a `character` cleaned of all modifications.
#' @noRd
#' @examples
#' .proforma_clean_sequences(
#'    c("EM[+15.9949]EVEES[+79.9663]PEK",
#'      "EM[+15.995]EVEES[-18.01]PEK")
#' )
.proforma_clean_sequences <- function(x) {
    gsub(pattern = "\\[[^]]*\\]|<[^>]*>", "", x)
}

#' Extract delta masses
#'
#' @name proforma-parser
#' @param x `character`, ProForma sequence.
#' @return `list`, a `list` of `doubles` representing the delta masses for each
#' sequence.
#' @noRd
#' @examples
#' .proforma_delta_masses(
#'    c("EM[+15.9949]EVEES[+79.9663]PEK",
#'      "EM[+15.995]EVEES[-18.01]PEK")
#' )
.proforma_delta_masses <- function(x) {
    rx <- gregexpr(
        pattern =
            "(?<AA>[A-Z])(?:\\[[GMURX]?:?)?(?<DeltaMass>[+-][0-9.]+)?(?:\\])?",
        text = x,
        perl = TRUE
    )

    mapply(function(sequence, start, matched_length) {
        mod <- as.double(
            substring(sequence, start, start + matched_length - 1L)
        )
        mod[is.na(mod)] <- 0
        mod
    },
        sequence = x,
        start =
            lapply(rx, function(r)attr(r, "capture.start")[, "DeltaMass"]),
        matched_length =
            lapply(rx, function(r)attr(r, "capture.length")[, "DeltaMass"]),
        SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
}

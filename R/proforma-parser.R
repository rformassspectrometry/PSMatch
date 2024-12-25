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
        pattern = "(?<=\\[)[GMURX]?:?[+-][0-9.]+(?=\\])",
        text = x,
        perl = TRUE
    )
    mapply(function(sequence, start, matched_length, n) {
        if (any(matched_length < 0))
            return(double(n))

        # add 2 for the surrounding "[" and "]"
        matched_length2 <- matched_length + 2L
        n_clean <- n - sum(matched_length2, na.rm = TRUE)
        masses <- double(n_clean)

        # subtract 2 for the "[" and the previous amino acid position
        masses[
            (start -
             cumsum(c(2L, matched_length2[-length(matched_length2)]) ))
        ] <- as.double(
            gsub(
                "^[GMURX]:",
                "",
                substring(sequence, start, start + matched_length - 1L)
            )
        )
        masses
    },
        sequence = x,
        start = rx,
        matched_length = lapply(rx, attr, "match.length"),
        n = nchar(x),
        SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
}

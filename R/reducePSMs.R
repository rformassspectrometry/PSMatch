##' @param x An instance of class `PSM`.
##'
##' @param k A `vector` or `factor` of length equal to `nrow(x)` that
##'     defines the primary key used to reduce `x`. This typically
##'     corresponds to the spectrum identifier, often names
##'     `spectrumID`.
##'
##' @return `reducePSMs()` returns a reduced version of the `x` input.
##'
##' @export reducePSMs
##'
##' @name PSM
reducePSMs <- function(x, k) {
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

##' @param object An instance of class `PSM`.
##'
##' @param k A `vector` or `factor` of length equal to `nrow(x)` that
##'     defines the primary key used to reduce `x`. This typically
##'     corresponds to the spectrum identifier. The defauls is to use
##'     the spectrum PSM variable.
##'
##' @return `reducePSMs()` returns a reduced version of the `x` input.
##'
##' @export reducePSMs
##'
##' @name PSM
reducePSMs <- function(object,
                       k = object[[psmVariables(object)["spectrum"]]]) {
    object <- QFeatures::reduceDataFrame(object, k)
    n <- ncol(object)
    for (i in seq_along(object)) {
        .x <- object[[i]]
        class_x <- class(.x)
        .x <- sapply(.x, unique)
        if (is.list(.x))
            .x <- as(.x, class_x)
        else .x <- unname(.x)
        object[[i]] <- .x
    }
    metadata(object)[["reduced"]] <- TRUE
    as(object, "PSM")
}

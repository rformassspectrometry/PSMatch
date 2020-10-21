##' Convert a `vector` of characters to camel case by replacing dots
##' by captial letters.
##'
##' @title Convert to camel case by replacing dots by captial letters
##' 
##' @param x A `character` to be transformed to camel case.
##' 
##' @param prefix An optional `character` of length one. Any
##'     additional elements are ignores.
##' 
##' @return A `character` of same length as `x`.
##' 
##' @author Laurent Gatto
##'
##' @noRd
##' 
##' @examples
##' nms <- c("aa.foo", "ab.bar")
##' psm:::makeCamelCase(nms)
##' psm:::makeCamelCase(nms, prefix = "x")
makeCamelCase <- function(x, prefix) {
    if (!missing(prefix))
        x <- paste(prefix[1], x, sep = ".")
    gsub('\\.(\\w?)', '\\U\\1', x, perl = TRUE)
}



##' This function produces the opposite as the `stringsAsFactors`
##' argument in the `data.frame` or `read.table` functions; it
##' converts `factors` columns to `characters`.
##'
##' @title Converts factors to strings
##' 
##' @param x A `data.frame`
##' 
##' @return A `data.frame` where `factors` are converted to
##'     `characters`.
##' 
##' @author Laurent Gatto
##'
##' @noRd
##' 
##' @examples
##' data(iris)
##' str(iris)
##' str(psm:::factorsAsStrings(iris))
factorsAsStrings <- function(x) {
    x <- lapply(x,
                   function(xx) {
                       if (is.factor(xx)) as.character(xx)
                       else xx
                   })
    data.frame(x, stringsAsFactors = FALSE)
}

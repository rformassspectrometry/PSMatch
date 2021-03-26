##' @exportClass PSM
setClass("PSM",
         contains = "DFrame")

setMethod("show", "PSM",
          function(object) {
              cl <- classNameForDisplay(object)
              if (!is.null(metadata(object)[["reduced"]]) &&
                  metadata(object)[["reduced"]])
                  cl <- paste("Reduced", cl)
              cat(cl, "with", nrow(object), "rows and",
                  ncol(object), "columns.\n")
              if (ncol(object) <= 4)
                  cat("names(", ncol(object), "): ",
                      paste(names(object), collapse = " "), "\n",
                      sep = "")
              else
                  cat("names(", ncol(object), "): ",
                      paste(names(object)[1:2], collapse = " "),
                      " ... ", paste(names(object)[(ncol(object)-1):ncol(object)],
                                   collapse = " "), "\n",
                      sep = "")
              invisible(NULL)
          })

##' @export
##'
##' @importFrom ProtGenerics adjacencyMatrix
##'
##' @rdname PSM
setMethod("adjacencyMatrix", "PSM",
          function(object) {
              if (is.na(psmVariables(object)[["protein"]]) |
                  is.na(psmVariables(object)[["peptide"]]))
                  stop("Please define the 'protein' and 'peptide' PSM variables.")
              if (!psmVariables(object)[["protein"]] %in% names(object) |
                  !psmVariables(object)[["peptide"]] %in% names(object))
                  stop("PSM variables 'protein' and 'peptide' must be defined.")
              vec <- object[[psmVariables(object)[["protein"]]]]
              names(vec) <- object[[psmVariables(object)[["peptide"]]]]
              ## NB: Note that we ignore any score here and always
              ## return a binary matrix. Use makeAdjacencyMatrix() for
              ## these features.
              makeAdjacencyMatrix(vec, binary = TRUE)
          })


##' @export
##'
##' @importFrom ProtGenerics adjacencyMatrix
##'
##' @rdname ConnectedComponents
setMethod("adjacencyMatrix", "ConnectedComponents",
          function(object) object@adjMatrix)

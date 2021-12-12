##' @export
##'
##' @importFrom ProtGenerics adjacencyMatrix
##'
##' @rdname PSM
setMethod("adjacencyMatrix", "PSM",
          function(object) {
              if (is.na(psmVariables(object)[["protein"]]) | is.na(psmVariables(object)[["peptide"]]))
                  stop("Please define the 'protein' and 'peptide' PSM variables.")
              vec <- object[[psmVariables(object)[["protein"]]]]
              names(vec) <- object[[psmVariables(object)[["peptide"]]]]
              makeAdjacencyMatrix(vec)
          })


##' @export
##'
##' @importFrom ProtGenerics adjacencyMatrix
##'
##' @rdname ConnectedComponents
setMethod("adjacencyMatrix", "ConnectedComponents",
          function(object) object@adjMatrix)

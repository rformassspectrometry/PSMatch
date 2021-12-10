setClass("ConnectedComponents",
         slots = c(adjMatrix = "Matrix",
                   ccMatrix = "Matrix",
                   ccList = "List",
                   ccPeptides = "List"))

ConnectedComponents <- function(adj) {
    getCCpeptides <- function(cc, adj) {
        res <- adj[, cc, drop = FALSE]
        res <- res[rowSums(res) > 0, , drop = FALSE]
        res
    }
    cc <- tcrossprod(t(adj))
    n <- ncol(cc)
    cc_pep <- cc_list <- vector("list", length = n)
    for (i in seq_len(n-1)) {
        j <- i:n
        k <- which(cc[i, j] != 0)
        cc_list[[i]] <- colnames(adj)[j][k]
        cc_pep[[i]] <- rownames(getCCpeptides(cc_list[[i]], adj))
    }
    sel <- lengths(cc_list) > 1
    new("ConnectedComponents",
        adjMatrix = adj,
        ccMatrix = cc,
        ccList = List(cc_list[sel]),
        ccPeptides = List(cc_pep[sel]))
}


setMethod("show", "ConnectedComponents",
          function(object) {
              cat(sprintf("An instance of class %s", class(object)), "\n")
              cat(" Number of componenents ", nrow(object@ccMatrix), "\n")
              cat(" CC with size > 1:\n")
              tab <- table(lengths(object@ccList))
              msg <- strwrap(paste(paste0(tab, "(", names(tab), ")"),
                                   collapse = " "))
              message(paste(" ", msg, collapse = "\n"))
          })

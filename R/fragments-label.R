##' @name labelFragments
##'
##' @title labels MS2 Fragments
##'
##' @description
##' Creates a list of annotations based on `calculateFragments` results.
##'
##' @param x An instance of class `Spectra` of length 1, containing a
##'     spectra variable `"sequence"` with a `character(1)`
##'     representing a valid peptide sequence.
##'
##' @param ppm m/z relative acceptable difference (in ppm) for peaks to be
##'     considered matching (see [MsCoreUtils::closest()] for more details).
##'
##' @param tolerance absolute acceptable difference of m/z values for peaks to
##'     be considered matching (see [MsCoreUtils::closest()] for more details).
##'
##' @param what `character(1)`, one of `"ion"` (default) or `"mz"`, defining
##'     whether labels should be fragment ions, , or their m/z values. If the
##'     latter, then the m/z values are named with the ion labels.
##'
##' @param ... additional parameters (except `verbose`) passed to
##'     [calculateFragments()] to calculate fragment m/z values to be
##'     added to the spectra in `x`.
##'
##' @return Return a `list()` of `character()` with fragment ion labels. The
##' elements are named after the peptide they belong to (variable
##' modifications included).
##'
##' @importFrom MsCoreUtils closest
##'
##' @export
##'
##' @author Johannes Rainer, Guillaume Deflandre, Sebastian Gibb, Laurent Gatto
##'
##' @examples
##' library("Spectra")
##'
##' ## Load the PSM and spectra data
##' data("psmBoekweg")
##' data("spBoekweg")
##'
##' ## The spectra need a 'sequence' variable for `calculateFragments()`
##' ## Make sure you can join both psms and spectra using `joinSpectraData`:
##' head(psmBoekweg$pkey <- paste0(basename(psmBoekweg$filename),
##'                                sub("^.+scan=", "::", psmBoekweg$scannr)))
##' head(spBoekweg$pkey <- paste0(basename(spBoekweg$dataOrigin),
##'                               sub("^.+scan=", "::", spBoekweg$spectrumId)))
##'
##' sp <- Spectra::joinSpectraData(spBoekweg, psmBoekweg, by.x = "pkey")
##' sp <- sp[sp$msLevel == 2]
##'
##' ## Add carbamidomethylation or other modifications if need be
##' ## See ?PTMods::addFixedModifications
##' (seq <- psmVariables(psmBoekweg)[["peptide"]])
##' (fdr <- psmVariables(psmBoekweg)[["fdr"]])
##'
##' head(sp[["sequence"]] <- addFixedModifications(sp[[seq]]))
##'
##' ## The fragment ion labels
##' labelFragments(sp[1])
##'
##' ## The fragment mz labels
##' labelFragments(sp[1], what = "mz")
##'
##' ## Pass additional parameters to calculateFragments using a PTMods modified sequence
##' sp_mod <- sp[1]
##' sp_mod$sequence <- PTMods::addFixedModifications(sp_mod$sequence,
##'                                                   fixedModifications = c(Nterm = 5))
##' labelFragments(sp_mod, type = c("a", "b", "x", "y"))
##'
##' ## Annotate the spectum with the fragment labels
##' plotSpectra(sp[1], labels = labelFragments, labelPos = 3)
##'
##' ## By default used in `plotSpectraPTM()`.
##' plotSpectraPTM(sp[1])
labelFragments <- function(x, tolerance = 0, ppm = 20,
                           what = c("ion", "mz"), ...) {
    stopifnot(requireNamespace("Spectra"))
    stopifnot(inherits(x, "Spectra"))
    what <- match.arg(what)
    super_labels <- vector("list", length = length(x))
    k <- integer()
    v <- peaksData(x)

    for (j in seq_along(x)) {
        stopifnot("sequence" %in% Spectra::spectraVariables(x[j]))
        y <- Spectra::spectraData(x[j])[["sequence"]]
        x_data <- v[[j]]
        y_data <- calculateFragments(y, verbose = FALSE, ...)

        y_data <- split(y_data, y_data$peptide)

        labels <- vector("list", length = length(y_data))
        names(labels) <- names(y_data)

        for (i in seq_along(y_data)) {
            k <- c(k, j)
            y_data[[i]] <- y_data[[i]][order(y_data[[i]]$mz), ]
            idy_all <- MsCoreUtils::closest(x_data[, "mz"],
                y_data[[i]][, "mz"],
                tolerance = tolerance,
                ppm = ppm
            )
            idx <- which(!is.na(idy_all))
            idy <- idy_all[idx]
            if (what == "ion") {
                labels[[i]] <- rep(NA_character_, nrow(x_data))
                labels[[i]][idx] <- y_data[[i]][idy, "ion"]
            } else { ## mz
                labels[[i]] <- rep(NA_real_, nrow(x_data))
                labels[[i]][idx] <- y_data[[i]][idy, "mz"]
            }
        }
        super_labels[[j]] <- labels
    }
    super_labels <- unlist(super_labels, recursive = FALSE)
    attr(super_labels, "group") <- k
    super_labels
}

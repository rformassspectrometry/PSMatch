##' @title Adds MS2 Fragments
##'
##' @param x An instance of class `Spectra` of length 1, containing a
##'     spectra variable `"sequence"` with a `character(1)`
##'     representing a valid peptide sequence.
##'
##' @param ppm m/z relative acceptable difference (in ppm) for peaks
##'     to be considered matching (see [MsCoreUtils::common()] for
##'     more #' details).
##'
##' @param tolerance absolute acceptable difference of m/z values for
##'     peaks to be considered matching (see [MsCoreUtils::common()]
##'     for more details).
##'
##' @param ... additional parameters (except `verbose`) passed to
##'     [calculateFragments()] to calculate fragment m/z values to be
##'     added to the spectra in `x`.
##'
##' @return Return a `list()` of `character()` with fragment ion labels. The 
##' elements are named after the peptide they belong to (modifications included).
##'
##' @importFrom MsCoreUtils common
##'
##' @export
##'
##' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
##'
##' @examples
##'
##' library("Spectra")
##'
##' sp <- DataFrame(msLevel = 2L, rtime = 2345, sequence = "SIGFEGDSIGR")
##' sp$mz <- list(c(100.048614501953, 110.069030761719, 112.085876464844,
##'                 117.112571716309, 158.089569091797, 163.114898681641,
##'                 175.117172241211, 177.098587036133, 214.127075195312,
##'                 232.137542724609, 233.140335083008, 259.938415527344,
##'                 260.084167480469, 277.111572265625, 282.680786132812,
##'                 284.079437255859, 291.208282470703, 315.422576904297,
##'                 317.22509765625, 327.2060546875, 362.211944580078,
##'                 402.235290527344, 433.255004882812, 529.265991210938,
##'                 549.305236816406, 593.217041015625, 594.595092773438,
##'                 609.848327636719, 631.819702148438, 632.324035644531,
##'                 632.804931640625, 640.8193359375, 641.309936523438,
##'                 641.82568359375, 678.357238769531, 679.346252441406,
##'                 688.291259765625, 735.358947753906, 851.384033203125,
##'                 880.414001464844, 881.40185546875, 919.406433105469,
##'                 938.445861816406, 1022.56658935547, 1050.50415039062,
##'                 1059.82800292969, 1107.52734375, 1138.521484375,
##'                 1147.51538085938, 1226.056640625))
##' sp$intensity <- list(c(83143.03, 65473.8, 192735.53, 3649178.5,
##'                        379537.81, 89117.58, 922802.69, 61190.44,
##'                        281353.22, 2984798.75, 111935.03, 42512.57,
##'                        117443.59, 60773.67, 39108.15, 55350.43,
##'                        209952.97, 37001.18, 439515.53, 139584.47,
##'                        46842.71, 1015457.44, 419382.31, 63378.77,
##'                        444406.66, 58426.91, 46007.71, 58711.72,
##'                        80675.59, 312799.97, 134451.72, 151969.72,
##'                        3215457.75, 1961975, 395735.62, 71002.98,
##'                        69405.73, 136619.47, 166158.69, 682329.75,
##'                        239964.69, 242025.44, 1338597.62, 50118.02,
##'                        1708093.12, 43119.03, 97048.02, 2668231.75,
##'                        83310.2, 40705.72))
##' sp <- Spectra(sp)
##'
##' ## The fragment ion labels
##' addFragments(sp)
##'
##' ## Annotate the spectum with the fragment labels
##' plotSpectra(sp, labels = addFragments, labelPos = 3)
addFragments <- function(x, tolerance = 0, ppm = 20, ...) {
    stopifnot(requireNamespace("Spectra"))
    stopifnot(inherits(x, "Spectra"))
    super_labels <- vector("list", length = length(x))
    
    for (j in seq_along(x)) {
        stopifnot("sequence" %in% Spectra::spectraVariables(x[j]))
        y <- Spectra::spectraData(x[j])[["sequence"]]
        x_data <- Spectra::peaksData(x[j])[[1L]]
        y_data <- calculateFragments(y, verbose = FALSE, ...)
        
        y_data <- split(y_data, y_data$peptide)
        
        labels <- vector("list", length = length(y_data))
        names(labels) <- names(y_data)
        
        for (i in seq_along(y_data)) {
            y_data[[i]] <- y_data[[i]][order(y_data[[i]]$mz), ]
            idx <- which(MsCoreUtils::common(x_data[, "mz"], 
                                             y_data[[i]][, "mz"],
                                             tolerance = tolerance,
                                             ppm = ppm))
            idy <- which(MsCoreUtils::common(y_data[[i]][, "mz"], 
                                             x_data[, "mz"],
                                             tolerance = tolerance, 
                                             ppm = ppm))
            
            labels[[i]] <- rep(NA_character_, nrow(x_data))
            labels[[i]][idx] <- y_data[[i]][idy, "ion"]
        }
        super_labels[[j]] <- labels
    }
    unlist(super_labels, recursive = FALSE)
}
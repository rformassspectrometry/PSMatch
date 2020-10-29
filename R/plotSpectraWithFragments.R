##' Ths
##'
##' @title Plots an MS2 Spectrum with Fragments
##' 
##' @param x An instance of class `Spectra`.
##' @param y `character(1)` with a peptides sequence.
##' @param xlab character(1)` with the label for the x-axis (by
#'     default `xlab = "m/z"`).
##' @param ylab `character(1)` with the label for the y-axis (by
#'     default `ylab = "intensity"`).
##' @param type `character(1)` specifying the type of plot. See
#'     [plot.default()] for details. Defaults to `type = "h"` which
#'     draws each peak as a line.
##' @param xlim `numeric(2)` defining the x-axis limits. The range of
#'     m/z values are used by default.
##' @param ylim `numeric(2)` defining the y-axis limits. The range of
#'     intensity values are used by default.
##' @param main `character(1)` with the title for the plot. By default
#'     the peptide sequence `y` is used.
##' @param col color to be used to draw the peaks. Should be either of
#'     length 1, or equal to the number of spectra (to plot each
#'     spectrum in a different color) or be a `list` with colors for
#'     each individual peak in each spectrum.
##' @param labelCex `numeric(1)` giving the amount by which the text
#'     should be magnified relative to the default. See parameter
#'     `cex` in [par()].
##' @param labelSrt `numeric(1)` defining the rotation of the
#'     label. See parameter `srt` in [text()].
##' @param labelAdj see parameter `adj` in [text()].
##' @param labelPos see parameter `pos` in [text()].
##' @param labelOffset see parameter `offset` in [text()].
##' @param axes `logical(1)` whether (x and y) axes should be drawn.
##' @param frame.plot `logical(1)` whether a box should be drawn
#'     around the plotting area.
##' @param ppm m/z relative acceptable difference (in ppm) for peaks
#'     to be considered matching (see [MsCoreUtils::common()] for more
#'     details).
##' @param tolerance absolute acceptable difference of m/z values for
#'     peaks to be considered matching (see [MsCoreUtils::common()]
#'     for more details).
##' @param matchCol color for matching peaks.
##' @param matchLwd line width (`lwd`) to draw matching peaks. See
#'     [par()] for more details.
##' @param matchLty line type (`lty`) to draw matching peaks. See
#'     [par()] for more details.
##' @param matchPch point character (`pch`) to label matching
##'     peaks. Defaults to `matchPch = 16`, set to `matchPch = NA` to
##'     disable. See [par()] for more details
##' @param ... additional parameters to be passed to [plot.default()].
##' 
##' @return Used to create a plot.
##'
##' @importFrom graphics plot.xy
##'
##' @importFrom grDevices xy.coords
##'
##' @export 
##' 
##' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
##'
##' @examples
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
##' sp <- Spectra::Spectra(sp)
##' plotSpectraWithFragments(sp, sp$sequence, labelPos = 3)
plotSpectraWithFragments <- function(x, y, xlab = "m/z", ylab = "intensity",
                                    type = "h", xlim = numeric(),
                                    ylim = numeric(),
                                    main = y,
                                    col = "#00000080", 
                                    labelCex = 1, labelSrt = 0,
                                    labelAdj = NULL, labelPos = NULL,
                                    labelOffset = 0.5, axes = TRUE,
                                    frame.plot = axes, ppm = 20, tolerance = 0,
                                    matchCol = "#80B1D3", matchLwd = 1,
                                    matchLty = 1, matchPch = 16, ...) {
    stopifnot(requireNamespace("Spectra"))
    stopifnot(requireNamespace("MsCoreUtils"))
    isValidSequence <- !missing(y) && !is.na(y) && nchar(y)
    stopifnot(isValidSequence)
    
    ## Prepare x and y data
    x_data <- Spectra::peaksData(x)[[1L]]
    y_data <- calculateFragments(y, verbose = FALSE)
    y_data <- y_data[order(y_data$mz), ]

    ## Find common peaks and prepare annotations
    idx <- which(MsCoreUtils::common(x_data[, "mz"], y_data[, "mz"],
                                     tolerance = tolerance, ppm = ppm))
    idy <- which(MsCoreUtils::common(y_data[, "mz"], x_data[, "mz"],
                                     tolerance = tolerance, ppm = ppm))
    labels <- rep(NA_character_, nrow(x_data))
    labels[idx] <- y_data[idy, "ion"]

    ## Plot spectrum with fragment annotation
    Spectra:::.plot_single_spectrum(x, add = FALSE, type = type,
                                    xlim = xlim, ylim = ylim, main = main,
                                    xlab = xlab, ylab = ylab,
                                    col = col[[1L]], labels = labels,
                                    labelCex = labelCex, labelCol = matchCol,
                                    labelSrt = labelSrt, labelAdj = labelAdj,
                                    labelPos = labelPos,
                                    labelOffset = labelOffset, ...)
    
    ## Highlight common peaks
    if (length(idx)) {
        plot.xy(xy.coords(x_data[idx, "mz"], x_data[idx, "intensity"]),
                type = "h", col = matchCol, lwd = matchLwd, ...)
        plot.xy(xy.coords(x_data[idx, "mz"], x_data[idx, "intensity"]),
                type = "p", col = matchCol, pch = matchPch, ...)
    }    
}

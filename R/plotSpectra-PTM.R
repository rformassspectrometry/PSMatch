##' @title Function to plot MS/MS spectra with PTMs
##'
##' @description
##'
##' `plotSpectraPTM()` creates annotated visualisations of MS/MS spectra,
##' designed to explore fragment identifications and post-translational
##' modifications (PTMs).
##'
##' `plotSpectraPTM()` plots a spectrum's m/z values on the x-axis and
##' corresponding intensities on the y-axis, labeling the peaks according to
##' theoretical fragment ions (e.g., b, y, a, c, x, z) computed using
##' `labelFragments()` and `calculateFragments()`.
##'
##' @param x a `Spectra()` object.
##'
##' @param deltaMz `logical(1L)` If `TRUE`, adds an additional plot showing the
##'     difference of mass over charge between matched oberved and theoretical
##'     fragments in parts per million. Does not yet support modifications. The
##'     matching is based on `calculateFragments()` and needs a 'sequence'
##'     variable in `spectraVariables(x)`. Default is set to `TRUE`.
##'
##' @param ppm `integer(1L)` Sets the limits of the delta m/z plot and is passed
##'     to `labelFragments()`.
##'
##' @param xlab `character(1)` with the label for the x-axis (by default `xlab =
##'     "m/z"`).
##'
##' @param ylab `character(1)` with the label for the y-axis (by default `ylab =
##'     "intensity"`).
##'
##' @param xlim `numeric(2)` defining the x-axis limits. The range of m/z values
##'     are used by default.
##'
##' @param ylim `numeric(2)` defining the y-axis limits. The range of intensity
##'     values are used by default.
##'
##' @param main `character(1)` with the title for the plot. By default the
##'     spectrum's MS level and retention time (in seconds) is used.
##'
##' @param col Named `character(4L)`. Colors for the labels, the character names
##'     need to be "b", "y", "acxz" and "other", respectively for the b-ions,
##'     y-ions, a,c,x,z-ions and the unidentified fragments.
##'
##' @param labelCex `numeric(1)` giving the amount by which the text should be
##'     magnified relative to the default. See parameter `cex` in `par()`.
##'
##' @param labelSrt `numeric(1)` defining the rotation of the label.  See
##'     parameter `srt` in `text()`.
##'
##' @param labelAdj see parameter `adj` in `text()`.
##'
##' @param labelPos see parameter `pos` in `text()`.
##'
##' @param labelOffset see parameter `offset` in `text()`.
##'
##' @param asp for `plotSpectra()`, the target ratio (columns / rows) when
##'     plotting mutliple spectra (e.g. for 20 spectra use asp = 4/5 for 4
##'     columns and 5 rows or asp = 5/4 for 5 columns and 4 rows; see
##'     `grDevices::n2mfrow()` for details).  If `deltaMz` is `TRUE`, `asp` is
##'     ignored.
##'
##' @param minorTicks `logical(1L)`. If `TRUE`, minor ticks are added to the
##'     plots.  Default is set to `TRUE`.
##'
##' @param ... additional parameters to be passed to the `labelFragments()`
##'     function.
##'
##' @importFrom graphics layout par
##'
##' @importFrom grDevices n2mfrow
##'
##' @importFrom Spectra spectraVariables
##'
##' @author Johannes Rainer, Sebastian Gibb, Guillaume Deflandre, Laurent Gatto
##'
##' @seealso [Spectra::plotSpectra()]
##'
##' @return Creates a plot depicting an MS/MS-MS spectrum.
##'
##' @export
##'
##' @examples
##' library("Spectra")
##'
##' sp <- DataFrame(msLevel = 2L, rtime = 2345, sequence = "SIGFEGDSIGR")
##' sp$mz <- list(c(75.048614501953, 81.069030761719, 86.085876464844,
##'                 88.039, 158.089569091797, 163.114898681641,
##'                 173.128, 177.098587036133, 214.127075195312,
##'                 232.137542724609, 233.140335083008, 259.938415527344,
##'                 260.084167480469, 277.111572265625, 282.680786132812,
##'                 284.079437255859, 291.208282470703, 315.422576904297,
##'                 317.22509765625, 327.2060546875, 362.211944580078,
##'                 402.235290527344, 433.255004882812, 534.258783,
##'                 549.305236816406, 593.217041015625, 594.595092773438,
##'                 609.848327636719, 631.819702148438, 632.324035644531,
##'                 632.804931640625, 640.8193359375, 641.309936523438,
##'                 641.82568359375, 678.357238769531, 679.346252441406,
##'                 706.309623, 735.358947753906, 851.384033203125,
##'                 880.414001464844, 881.40185546875, 906.396433105469,
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
##'                        1961975, 69405.76, 395735.62, 71002.98,
##'                        3215457.75, 136619.47, 166158.69, 682329.75,
##'                        239964.69, 242025.44, 1338597.62, 50118.02,
##'                        1708093.12, 43119.03, 97048.02, 2668231.75,
##'                        83310.2, 40705.72))
##' sp <- Spectra(sp)
##'
##' ## Annotate the spectum with the fragment labels
##' plotSpectraPTM(sp, main = "An example of a annotated plot")
##'
##' ## Annotate the spectrum without the delta m/z plot
##' plotSpectraPTM(sp, deltaMz = FALSE)
##'
##' ## Annotate the spectrum with different ion types
##' plotSpectraPTM(sp, type = c("a", "b", "x", "y"))
##'
##' ## Annotate the spectrum with variable modifications
##' plotSpectraPTM(sp, variable_modifications = c(R = 49.469))
##'
##' ## Annotate multiple spectra at a time
##' plotSpectraPTM(c(sp,sp), variable_modifications = c(R = 490469))
##'
##' ## Color the peaks with different colors
##' plotSpectraPTM(sp, col = c(y = "red", b = "blue", acxy = "chartreuse3", other = "black"))
plotSpectraPTM <- function(x, deltaMz = TRUE, ppm = 20,
                           xlab = "m/z", ylab = "intensity [%]",
                           xlim = numeric(), ylim = numeric(),
                           main = character(),
                           col = c(y = "darkred",
                                   b = "darkblue",
                                   acxy = "darkgreen",
                                   other = "grey40"),
                           labelCex = 1, labelSrt = 0,
                           labelAdj = NULL, labelPos = 3, labelOffset = 0.5,
                           asp = 1, minorTicks = TRUE,
                           ...) {
    if (!("sequence" %in% Spectra::spectraVariables(x))) {
        stop("Missing 'sequence' in Spectra::spectraVariables(x)")
    }
    nsp <- length(x)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    if (length(main) != nsp) main <- rep(main[1], nsp)

    labels <- labelFragments(x, ppm = ppm, what = "ion", ...)

    if (deltaMz) { ## Generate deltaMzData labels for .plot_single_spectrum_PTM
        deltaMzData <- labelFragments(x, ppm = ppm, what = "mz", ...)
        layout_matrix <- .make_layout_matrix(length(labels))
        layout(layout_matrix,
               heights = rep(c(5, 1), length.out = nrow(layout_matrix)))
    } else {
      par(mfrow = n2mfrow(length(labels), asp = asp))
      deltaMzData <- NULL
    }
    spectrum_number <- attr(labels, "group")

    for (i in seq_along(spectrum_number)) {
        .plot_single_spectrum_PTM(x[spectrum_number[i]], xlab = xlab,
                                  ylab = ylab, xlim = xlim,
                                  ylim = ylim,
                                  main = main[spectrum_number[i]],
                                  col = col, labels = labels[i],
                                  labelCex = labelCex, labelSrt = labelSrt,
                                  labelAdj = labelAdj, labelPos = labelPos,
                                  labelOffset = labelOffset,
                                  minorTicks = minorTicks,
                                  deltaMzData = deltaMzData[[i]],
                                  ppm = ppm, old_par = old_par)
    }
}


##' @description
##'
##' Plot a single spectrum (m/z on x against intensity on y) whilst labeling the
##' individual peaks based on `labelFragments()`.
##'
##' @author Johannes Rainer, Sebastian Gibb, Guillaume Deflandre, Laurent Gatto
##'
##' @importFrom graphics axis plot.new plot.window plot.xy strwidth
##'
##' @importFrom graphics text title par mtext abline points
##'
##' @importFrom grDevices xy.coords
##'
##' @importFrom Spectra peaksData spectraData
##'
##' @noRd
.plot_single_spectrum_PTM <- function(x, xlab = "m/z", ylab = "intensity",
                                      xlim = numeric(),
                                      ylim = numeric(), main = character(),
                                      col = c(y = "darkred",
                                              b = "darkblue",
                                              acxy = "darkgreen",
                                              other = "grey40"),
                                      labels = list(),
                                      labelCol = col, labelCex = 1, labelSrt = 0,
                                      labelAdj = NULL, labelPos = 3,
                                      labelOffset = 0.5, minorTicks = TRUE,
                                      ppm = 20, deltaMzData = NULL,
                                      old_par = old_par) {
    v <- peaksData(x)[[1L]]
    mzs <- v[, "mz"]
    ints <- v[, "intensity"]

    if (!length(xlim))
        suppressWarnings(xlim <- range(mzs, na.rm = TRUE))
    if (!length(ylim))
        suppressWarnings(
            ylim <- range(c(0, max(abs(ints), na.rm = TRUE))))
    if (any(is.infinite(xlim)))
        xlim <- c(0, 0)
    if (any(is.infinite(ylim)))
        ylim <- c(0, 0)
    if (!is.na(main)) {
      par(cex.main = 1.2, mar = c(4,4,1.5,2))
    } else par(mar = c(4,4,1,2))
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)

    peptide_sequence <- names(labels)
    labels <- labels[[1]]
    wdths <- max(strwidth(labels, cex = labelCex)) / 2
    usr_lim <- par("usr")
    ylim[2L] <- ylim[2L] + 0.6 * diff(ylim)
    xlim[1L] <- xlim[1L] - wdths
    xlim[2L] <- xlim[2L] + wdths
    plot.window(xlim = xlim, ylim = ylim)

    peakCol <- rep_len(col[["other"]], length(labels))
    peakCol[startsWith(labels, "b")] <- col[["b"]]
    peakCol[startsWith(labels, "y")] <- col[["y"]]
    peakCol[grepl("^[acxz]", labels)] <- col[["acxy"]]

    labelCol <- ifelse(grepl("b", labels), col[["b"]],
                ifelse(grepl("y", labels), col[["y"]], col[["acxy"]]))

    text(mzs, ints, labels = labels, adj = labelAdj, pos = labelPos,
         col = labelCol, cex = labelCex, srt = labelSrt, offset = labelOffset)

    .build_annotations(x, mzs, ints, col, labels, peptide_sequence)

    plot.xy(xy.coords(mzs, ints), type = "h", col = peakCol)

    major_ticks <- pretty(mzs, n = 8)
    axis(side = 1, lwd = 1, at = major_ticks, pos = 0,
         col.ticks = "grey45", col = "grey45")

    if (minorTicks) {
      nm <- length(major_ticks)
      ticks <- seq.int(
        major_ticks[1L], major_ticks[nm], length.out = 5L * (nm - 1L) + 1L)
      axis(side = 1, at = ticks[!ticks %in% major_ticks], labels = FALSE,
           tck = -0.01, col.ticks = "grey65", pos = 0)
    }

    axis_y_percent <- seq(0, 100, length.out = 5)

    axis(side = 2, las = 0, tck = -0.02,
         at = seq(0, max(abs(ints)), length.out = 5),
         labels = axis_y_percent)

    title(main = main, xlab = xlab, ylab = ylab)

    prefix <- "mzspec"
    run <- basename(spectraData(x)[["dataOrigin"]])
    scan <- paste0("scan: ", spectraData(x)[["scanIndex"]])
    rt <- paste0("rt: ", round(spectraData(x)[["rtime"]], 2))
    charge <- paste0("charge: ", spectraData(x)[["charge"]])
    seq_text <- paste0("sequence: ", peptide_sequence)

    mtext(paste(prefix, run, scan, rt, charge, seq_text, sep = "/"),
          cex = 0.85, line = -1)

    base_peak <- which.max(abs(ints))
    text(mzs[base_peak], ints[base_peak] * 0.60,
         paste0(formatC(ints[base_peak])),
         pos = 4, offset = 0.6, cex = 0.9, srt = 90)

    abline(h = 0, col = "grey45")

    if (!is.null(deltaMzData)) {
        deltaMzData <-
            ((mzs - deltaMzData) / deltaMzData) * 10^6
        par(mar = c(2, 4, 0, 2))

        true_hits <- !is.na(labels)

        plot(
            mzs[true_hits],
            deltaMzData[true_hits],
            col = labelCol[true_hits],
            type = "h", pch = 19,
            ann = FALSE, xaxt = "n",
            xlim = xlim, ylim = c(-ppm, ppm),
            lwd = 2)

        points(
            mzs[true_hits], deltaMzData[true_hits],
            col = labelCol[true_hits],
            type = "p", pch = 19, cex = 0.7)

        abline(h = 0, col = "#808080", lty = 2)
        title(ylab = "delta m/z\n[ppm]", cex.lab = 0.9, line = 2)
        axis(side = 1, lwd = 1, at = pretty(xlim, n = 8),
             col.ticks = "grey45", col = "grey45")
    }
}

##' Build layout matrix of plots to be printed
##'
##' @param n Number of spectra to be plotted
##'
##' @author Guillaume Deflandre
##'
##' @return A layout matrix
##'
##' @noRd
.make_layout_matrix <- function(n) {
    n_cols <- ceiling(sqrt(n))
    n_rows <- ceiling(n / n_cols)

    layout_matrix <- matrix(0, nrow = 2 * n_rows, ncol = n_cols)

    plot_index <- 1
    for (i in seq_len(n)) {
        row_block <- ((i - 1) %/% n_cols) * 2
        col_block <- ((i - 1) %% n_cols) + 1

        layout_matrix[row_block + 1, col_block] <- plot_index     # MSMS
        layout_matrix[row_block + 2, col_block] <- plot_index + 1 # delta m/z

        plot_index <- plot_index + 2
    }

    layout_matrix
}

##' Annotated sequence fragments split view
##'
##' @param x Spectra object
##'
##' @param peptide_sequence `character(1L)` The identified peptide sequence
##'
##' @param mzs mz values of peaks
##'
##' @param ints the intensity values of peaks
##'
##' @param labels labels associated to the peaks
##' 
##' @param col The colors of annotations
##'
##' @importFrom graphics segments
##'
##' @noRd
.build_annotations <- function(x,
                               mzs = mzs,
                               ints = ints,
                               col = col,
                               labels = labels,
                               peptide_sequence = peptide_sequence) {
  split_seq <- unlist(
    strsplit(
      peptide_sequence,
      "(?<=[A-Za-z])(?=[A-Z])|(?<=\\])(?=[A-Z])",
      perl = TRUE
    ))
  
  cleaned_peptide <- sapply(split_seq, function(amino_acid) {
    gsub("([A-Za-z])\\[[-+]?\\d+\\.?\\d*\\]", "[\\1]", amino_acid)
  })
  
  res_count <- length(cleaned_peptide)
  label_template <- rep(NA, res_count * 2 - 1)
  label_template[seq(1, length(label_template), by = 2)] <- cleaned_peptide
  
  peptide_annotation <- label_template
  
  .insert_ions <- function(prefix, idx) {
    label_template[seq(2, length(label_template), by = 2)] <- paste0(prefix, idx)
    return(label_template)
  }
  
  b_ion_list <- .insert_ions("b", seq_len(res_count - 1))
  y_ion_list <- .insert_ions("y", rev(seq_len(res_count - 1)))
  
  mz_center <- mean(range(mzs))
  mz_range <- diff(range(mzs)) * 0.4
  x_positions <- seq(mz_center - mz_range / 2, 
                     mz_center + mz_range / 2, 
                     length.out = length(label_template))
  
  ion_label_data <- data.frame(
    x = x_positions,
    residue = unlist(peptide_annotation),
    b_label = unlist(b_ion_list),
    y_label = unlist(y_ion_list),
    stringsAsFactors = FALSE
  )
  
  label_height <- max(abs(ints)) * 1.35
  label_spacing <- max(abs(ints)) / 15
  b_color <- col[["b"]]
  y_color <- col[["y"]]
  
  displayed_residues <- ion_label_data[!is.na(ion_label_data$residue), ]
  
  text(
    displayed_residues$x,
    label_height,
    displayed_residues$residue,
    cex = 1,
    adj = 0.5
  )
  
  b_ion_hits <- ion_label_data[ion_label_data$b_label %in% labels, ]
  if (nrow(b_ion_hits) > 0) {
    b_x <- b_ion_hits$x
    b_idx <- match(b_x, ion_label_data$x) - 1
    b_midpoints <- (b_x + ion_label_data$x[b_idx]) / 2
    b_ion_tags <- substring(b_ion_hits$b_label, 2)
    
    segments(
      b_x,
      label_height,
      b_x,
      label_height - label_spacing,
      col = b_color,
      lwd = 2
    )
    segments(
      b_midpoints,
      label_height - label_spacing * 1.5,
      b_x,
      label_height - label_spacing,
      col = b_color,
      lwd = 2
    )
    text(b_midpoints + (b_midpoints - b_x)/2,
         label_height - label_spacing * 1.5,
         b_ion_tags,
         adj = c(0.5, 1.3),
         col = b_color
    )}
  
  y_ion_hits <- ion_label_data[ion_label_data$y_label %in% labels, ]
  if (nrow(y_ion_hits) > 0) {
    y_x <- y_ion_hits$x
    y_idx <- match(y_x, ion_label_data$x) + 1
    y_midpoints <- (y_x + ion_label_data$x[y_idx]) / 2
    y_ion_tags <- substring(y_ion_hits$y_label, 2)
    
    segments(
      y_x,
      label_height,
      y_x,
      label_height + label_spacing,
      col = y_color,
      lwd = 2
    )
    segments(
      y_midpoints,
      label_height + label_spacing * 1.5,
      y_x,
      label_height + label_spacing,
      col = y_color,
      lwd = 2
    )
    text(y_midpoints + (y_midpoints - y_x)/2,
         label_height + label_spacing * 1.5,
         y_ion_tags,
         adj = c(0.5, -0.3),
         col = y_color
    )}
}





#' @title Function to plot MS/MS spectra with PTMs
#'
#' @description
#' Plot a single spectrum, annotate peaks, easily visualise sequence
#' fragmentations, ...
#' 
#' @name plotSpectraPTM
#'
#' @param x a `Spectra()` object.
#'
#' @param xlab `character(1)` with the label for the x-axis (by default `xlab = "m/z"`).
#'
#' @param ylab `character(1)` with the label for the y-axis (by default `ylab = "intensity"`).
#'
#' @param type `character(1)` specifying the type of plot. See `plot.default()`
#' for details. Defaults to `type = "h"` which draws each peak as a line.
#'
#' @param xlim `numeric(2)` defining the x-axis limits. The range of m/z values are used by default.
#'
#' @param ylim `numeric(2)` defining the y-axis limits. The range of intensity values are used by default.
#'
#' @param main `character(1)` with the title for the plot. By default the spectrum's
#'  MS level and retention time (in seconds) is used.
#'
#' @param col color to be used to draw the peaks. Should be either of length 1,
#' or equal to the number of spectra (to plot each spectrum in a different color)
#' or be a list with colors for each individual peak in each spectrum.
#'
#' @param labels allows to specify a label for each peak. Can be a character with
#'  length equal to the number of peaks, or, ideally, a function that uses one of
#'  the Spectra's variables (see examples below).
#'
#' @param labelCex `numeric(1)` giving the amount by which the text should be
#' magnified relative to the default. See parameter `cex` in `par()`.
#'
#' @param labelSrt `numeric(1)` defining the rotation of the label. See parameter `srt` in `text()`.
#'
#' @param labelAdj see parameter `adj` in `text()`.
#'
#' @param labelPos see parameter `pos` in `text()`.
#'
#' @param labelOffset see parameter `offset` in `text()`.
#'
#' @param labelCol color for the label(s).
#'
#' @param asp for `plotSpectra()`: the target ratio (columns / rows) when plotting
#' mutliple spectra (e.g. for 20 spectra use asp = 4/5 for 4 columns and 5 rows
#' or asp = 5/4 for 5 columns and 4 rows; see `grDevices::n2mfrow()` for details).
#'
#' @param ... additional parameters to be passed to the `plot.default()` function.
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#' 
#' @return Creates a plot depicting an MS/MS-MS spectrum
#'
#' @export
#' 
#' @examples
#' library("Spectra")
#' sp <- DataFrame(msLevel = 2L, rtime = 2345, sequence = "SIGFEGDSIGR")
#' sp$mz <- list(c(100.048614501953, 110.069030761719, 112.085876464844,
#'                 117.112571716309, 158.089569091797, 163.114898681641,
#'                 175.117172241211, 177.098587036133, 214.127075195312,
#'                 232.137542724609, 233.140335083008, 259.938415527344,
#'                 260.084167480469, 277.111572265625, 282.680786132812,
#'                 284.079437255859, 291.208282470703, 315.422576904297,
#'                 317.22509765625, 327.2060546875, 362.211944580078,
#'                 402.235290527344, 433.255004882812, 529.265991210938,
#'                 549.305236816406, 593.217041015625, 594.595092773438,
#'                 609.848327636719, 631.819702148438, 632.324035644531,
#'                 632.804931640625, 640.8193359375, 641.309936523438,
#'                 641.82568359375, 678.357238769531, 679.346252441406,
#'                 688.291259765625, 735.358947753906, 851.384033203125,
#'                 880.414001464844, 881.40185546875, 919.406433105469,
#'                 938.445861816406, 1022.56658935547, 1050.50415039062,
#'                 1059.82800292969, 1107.52734375, 1138.521484375,
#'                 1147.51538085938, 1226.056640625))
#' sp$intensity <- list(c(83143.03, 65473.8, 192735.53, 3649178.5,
#'                        379537.81, 89117.58, 922802.69, 61190.44,
#'                        281353.22, 2984798.75, 111935.03, 42512.57,
#'                        117443.59, 60773.67, 39108.15, 55350.43,
#'                        209952.97, 37001.18, 439515.53, 139584.47,
#'                        46842.71, 1015457.44, 419382.31, 63378.77,
#'                        444406.66, 58426.91, 46007.71, 58711.72,
#'                        80675.59, 312799.97, 134451.72, 151969.72,
#'                        3215457.75, 1961975, 395735.62, 71002.98,
#'                        69405.73, 136619.47, 166158.69, 682329.75,
#'                        239964.69, 242025.44, 1338597.62, 50118.02,
#'                        1708093.12, 43119.03, 97048.02, 2668231.75,
#'                        83310.2, 40705.72))
#' sp <- Spectra(sp)
#' ## Annotate the spectum with the fragment labels
#' plotSpectraPTM(sp, labels = addFragments2, labelPos = 3)
plotSpectraPTM <- function(x, xlab = "m/z", ylab = "intensity", type = "h",
                           xlim = numeric(), ylim = numeric(),
                           main = NULL, col = "#00000080",
                           labels = character(), labelCex = 1, labelSrt = 0,
                           labelAdj = NULL, labelPos = 3, labelOffset = 0.5,
                           labelCol = "#00000080", asp = 1, ...) {
    
    nsp <- length(x)
    old_par <- par(no.readonly = TRUE)
    if (nsp == 1)
        col <- list(col)
    if (length(col) != nsp)
        col <- rep(col[1], nsp)
    if (length(main) && length(main) != nsp)
        main <- rep(main[1], nsp)
    if (nsp > 1) {
        par(mfrow = n2mfrow(nsp, asp = asp))
    }
    
    ## only applicable for addFragments and not custom functions ...
    if (length(labels)) {
        if (is.function(labels)) {
            labels <- labels(x)
        }
        col <- rep(col[1], length(labels))
        par(mfrow = n2mfrow(length(labels), asp = asp))
        
        spectrum_number <- sapply(labels, function(label) attr(label, "spectrumNumber"))
        
        for (i in seq_along(spectrum_number)) {
            .plot_single_spectrum(x[spectrum_number[i]], xlab = xlab, ylab = ylab, 
                                  type = type, xlim = xlim, ylim = ylim, 
                                  main = main[spectrum_number[i]],
                                  col = col[[spectrum_number[i]]], 
                                  labels = labels[i],
                                  labelCex = labelCex, labelSrt = labelSrt,
                                  labelAdj = labelAdj, labelPos = labelPos,
                                  labelOffset = labelOffset, labelCol = labelCol,
                                  ...)
        }
    } else {
        for (i in seq_len(nsp))
            .plot_single_spectrum(x[i], xlab = xlab, ylab = ylab, type = type,
                                  xlim = xlim, ylim = ylim, main = main[i],
                                  col = col[[i]], labels = character(),
                                  labelCex = labelCex, labelSrt = labelSrt,
                                  labelAdj = labelAdj, labelPos = labelPos,
                                  labelOffset = labelOffset, labelCol = labelCol,
                                  ...)
    }
    on.exit(par(old_par))
    
}


#' @description
#'
#' Plot a single spectrum (m/z on x against intensity on y) with the optional
#' possibility to label the individual peaks.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @importFrom graphics axis box plot.new plot.window plot.xy strwidth
#'
#' @importFrom graphics text title
#'
#' @importFrom grDevices dev.flush dev.hold xy.coords
#'
#' @examples
#'
#' ints <- c(4.3412, 12, 8, 34, 23.4)
#' mzs <- c(13.453421, 43.433122, 46.6653553, 129.111212, 322.24432)
#'
#' df <- DataFrame(msLevel = 1L, rtime = 123.12)
#' df$mz <- list(mzs)
#' df$intensity <- list(ints)
#' sp <- Spectra(df)
#'
#' .plot_single_spectrum(sp, main = "hello")
#' .plot_single_spectrum(sp, bty = "n")
#' .plot_single_spectrum(sp, frame.plot = FALSE)
#'
#' .plot_single_spectrum(sp, col = 1:5)
#'
#' .plot_single_spectrum(sp, labels = 1:5, col = 1:5)
#'
#' .plot_single_spectrum(sp, labels = format(mz(sp)[[1]], digits = 5),
#'     labelPos = 2, labelOffset = 0.1, labelSrt = -30)
#' grid()
#' .plot_single_spectrum(sp, col = "red", type = "p", add = TRUE)
#'
#' .plot_single_spectrum(sp,
#'     labels = function(z) format(mz(z)[[1]], digits = 5),
#'     labelPos = 2, labelOffset = 0.1, labelSrt = -30)
#' grid()
#'
#' @noRd
.plot_single_spectrum <- function(x, xlab = "m/z", ylab = "intensity",
                                  type = "h", xlim = numeric(),
                                  ylim = numeric(),
                                  main = NULL,
                                  col = "#00000080", labels = character(),
                                  labelCol = col, labelCex = 1, labelSrt = 0,
                                  labelAdj = NULL, labelPos = 3,
                                  labelOffset = 0.5, add = FALSE,
                                  axes = TRUE, frame.plot = axes,
                                  orientation = 1, ...) {
    v <- peaksData(x)[[1L]]
    mzs <- v[, "mz"]
    ints <- orientation * v[, "intensity"]
    old_par <- par(no.readonly = TRUE)
    
    if (!length(xlim))
        suppressWarnings(xlim <- range(mzs, na.rm = TRUE))
    if (!length(ylim))
        suppressWarnings(
            ylim <- range(orientation * c(0, max(abs(ints), na.rm = TRUE))))
    if (any(is.infinite(xlim)))
        xlim <- c(0, 0)
    if (any(is.infinite(ylim)))
        ylim <- c(0, 0)
    if (!add) {
        if (length(main)) {
            par(mar = c(4,4,1.4,3))   
        } else {par(mar = c(4,4,0.5,3))}
        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
    }
    
    if (length(labels)) {
        peptide_sequence <- names(labels)
        labels <- labels[[1]]
        wdths <- max(strwidth(labels, cex = labelCex)) / 2
        usr_lim <- par("usr")
        if (orientation == 1) {
            ylim[2L] <- ylim[2L] + 0.7*diff(ylim)
        } else {
            ylim[1L] <- ylim[1L] - 0.7*diff(ylim)
        }
        xlim[1L] <- xlim[1L] - wdths
        xlim[2L] <- xlim[2L] + wdths
        if (!add)
            plot.window(xlim = xlim, ylim = ylim, ...)
        
        peakCol <- ifelse(grepl("b", labels), "blue",
                          ifelse(grepl("y", labels), "red", col))
        
        labelCol <- ifelse(grepl("b", labels), "darkblue",
                           ifelse(grepl("y", labels), "darkred", "darkgreen"))
        
        if (orientation == -1) labelPos = 1
        
        text(mzs, ints, labels = labels, adj = labelAdj, pos = labelPos,
             col = labelCol, cex = labelCex, srt = labelSrt,
             offset = labelOffset)
        
        .draw_psmanno(x, mzs, ints, labels, orientation, peptide_sequence)
        
        plot.xy(xy.coords(mzs, ints), type = type, col = peakCol, ...)
        
    }
    
    if (!add) {
        if (axes) {
            axis(side = 1, lwd = 0, ...)
            
            pretty_value <- pretty(seq(0, 1.07*max(abs(ints))*orientation, length.out = 8), n = 6)
            pretty_value_max <- pretty_value[length(pretty_value)]
            
            axis(side = 2, at = pretty_value, ...)
            
            axis_y_percent <- paste0(seq(0, 100, length.out = 5), "%")
            
            axis(side = 4, las = 0, tck = -0.02, 
                 at = seq(0, orientation*max(abs(ints)), length.out = 5),
                 labels = axis_y_percent,
                 col = "grey45", col.ticks = "grey45", col.axis = "grey45")
        }
        if (frame.plot)
            # box(...)
            title(main = main, xlab = xlab, ylab = ylab, ...)
    }
    plot.xy(xy.coords(mzs, ints), type = type, col = col, ...)
    
    prefix <- "mzspec"
    run <- basename(spectraData(x)[["dataOrigin"]])
    scan <- paste0("scan: ",spectraData(x)[["scanIndex"]])
    rt <- paste0("rt: ", round(spectraData(x)[["rtime"]], 2))
    charge <- paste0("charge: ", spectraData(x)[["charge"]])
    if (!length(labels)) {
        seq_text <- "sequence: unknown"
    } else {
        seq_text <- paste0("sequence: ", peptide_sequence)
    }
    
    text(min(mzs), max(abs(ints))*1.7*orientation,
         paste(prefix, run, scan, rt, charge, seq_text, sep = "/"),
         pos=4, offset=0, cex = 0.7)
    
    abline(h = 0, lwd = 0.5)
}


#' Annotated sequence fragments split view
#'
#' @param x Spectra object
#'
#' @noRd
#'
.draw_psmanno <- function(x, mzs = mzs, ints = ints,
                          labels = labels, orientation = orientation,
                          peptide_sequence = peptide_sequence){
    
    pep_seq <- unlist(strsplit(peptide_sequence,
                               "(?<=[A-Za-z])(?=[A-Z])|(?<=\\])(?=[A-Z])", perl = T))
    
    mod_pep_seq <- NULL
    for (i in seq_along(pep_seq)) {
        letter <- gsub("([A-Za-z])\\[[-+]?\\d+\\.?\\d*\\]", "[\\1]", pep_seq[i])
        mod_pep_seq <- c(mod_pep_seq, letter)
    }
    
    peptide <- c("-", mod_pep_seq, "-")
    
    peptide_list <- vector("list", length(peptide)*2-1) # set peptide + two dash
    peptide_list_b <- peptide_list # preset NULL for b ions
    peptide_list_y <- peptide_list # preset NULL for y ions
    
    # set peptide
    peptide_list[c(TRUE,FALSE)] <- as.list(peptide)
    peptide_list[c(FALSE,TRUE)] <- as.list(".")
    # set b ions
    peptide_list_b[c(FALSE,TRUE)] <-
        paste("b", seq(0,length(peptide)-2,by=1),sep="")
    # set y ions
    peptide_list_y[c(FALSE,TRUE)] <-
        paste("y", rev(seq(0,length(peptide)-2,by=1)), sep="")
    
    # set the width of AA sequence in the plot
    x_quarter<-seq(min(mzs),max(mzs), length.out = 20)[c(3,18)]  # seq from 5/20 to 17/20
    AA_pos <- seq(x_quarter[1], x_quarter[2], length.out=length(peptide_list))
    pos_start_dash <- x_quarter[1]  # x-axis value "-" at the front of peptides
    pos_end_dash   <- x_quarter[2] # x-axis value "-" at the end of peptides
    
    # integrate as data.table
    PSMlabel <-
        data.table::data.table(AA_pos = AA_pos, peptide = peptide_list,
                               bion = peptide_list_b, yion = peptide_list_y)
    
    peptide_height <- max(abs(ints))*1.4*orientation
    len_annoSpace <- max(abs(ints))/10*orientation
    
    b_ion_col <- "darkblue"
    y_ion_col <- "darkred"
    
    # Remove rows containing a dot (".") in any column
    PSMlabel_annots <-
        PSMlabel[!apply(PSMlabel, 1, function(row) any(row %in% c(".", "-"))), ]
    
    text(PSMlabel_annots$AA_pos, peptide_height,
         PSMlabel_annots$peptide, cex = 1, adj = 0.5 )
    
    # draw b ions between AA letters
    PSManno_bion <- subset(PSMlabel, PSMlabel$bion %in% labels)$AA_pos
    index <- match(PSManno_bion,PSMlabel$AA_pos)-1
    bion_xsmall <- (PSMlabel$AA_pos[index] + PSManno_bion)/2
    PSManno_bsmall <-
        substring(subset(PSMlabel, PSMlabel$bion %in% labels)$bion,2)
    
    # draw y ions between AA letters
    PSManno_yion <- subset(PSMlabel, PSMlabel$yion %in% labels)$AA_pos
    index <- match(PSManno_yion,PSMlabel$AA_pos) + 1
    yion_xlarge <- (PSMlabel$AA_pos[index] + PSManno_yion)/2
    PSManno_ysmall <-
        substring(subset(PSMlabel, PSMlabel$yion %in% labels)$yion, 2)
    
    if (orientation == 1) {
        if(length(PSManno_bion)>0){ # bion
            segments(PSManno_bion, rep(peptide_height, length(PSManno_bion)),
                     PSManno_bion,rep(peptide_height - len_annoSpace,length(PSManno_bion)),
                     col = b_ion_col, lwd = 2)
            segments(bion_xsmall, rep(peptide_height - len_annoSpace,
                                      length(PSManno_bion)), PSManno_bion,
                     rep(peptide_height - len_annoSpace, length(PSManno_bion)),
                     col = b_ion_col, lwd = 2)
            
            text((bion_xsmall + PSManno_bion)/2,
                 rep(peptide_height - len_annoSpace, length(PSManno_bion)),
                 PSManno_bsmall, cex = 1, adj = c(0.5, 1.3), col = b_ion_col)
        }
        if(length(PSManno_yion) > 0 ) {  # y ion
            segments(PSManno_yion, rep(peptide_height, length(PSManno_yion)),
                     PSManno_yion,rep(peptide_height + len_annoSpace,
                                      length(PSManno_yion)),
                     col = y_ion_col, lwd = 2)
            segments(yion_xlarge,
                     rep(peptide_height + len_annoSpace, length(PSManno_yion)),
                     PSManno_yion, rep(peptide_height + len_annoSpace,
                                       length(PSManno_yion)), 
                     col = y_ion_col, lwd = 2)
            
            text((yion_xlarge+PSManno_yion)/2,
                 rep(peptide_height + len_annoSpace, length(PSManno_yion)),
                 PSManno_ysmall, cex = 1, adj = c(0.5, -0.3),col = y_ion_col)
        }
    } else {
        if(length(PSManno_bion)>0){ # bion
            segments(PSManno_bion, rep(peptide_height, length(PSManno_bion)),
                     PSManno_bion,rep(peptide_height + len_annoSpace,length(PSManno_bion)),
                     col = b_ion_col, lwd = 2)
            segments(bion_xsmall, rep(peptide_height + len_annoSpace, length(PSManno_bion)), 
                     PSManno_bion, rep(peptide_height + len_annoSpace, length(PSManno_bion)),
                     col = b_ion_col, lwd = 2)
            
            text((bion_xsmall + PSManno_bion)/2,
                 rep(peptide_height + len_annoSpace, length(PSManno_bion)),
                 PSManno_bsmall, cex = 1, adj = c(0.5, 1.3), col = b_ion_col)
        }
        if(length(PSManno_yion) > 0 ) {  # y ion
            segments(PSManno_yion, rep(peptide_height, length(PSManno_yion)),
                     PSManno_yion, rep(peptide_height - len_annoSpace, length(PSManno_yion)),
                     col = y_ion_col, lwd = 2)
            segments(yion_xlarge, rep(peptide_height - len_annoSpace, length(PSManno_yion)),
                     PSManno_yion, rep(peptide_height - len_annoSpace, length(PSManno_yion)), 
                     col = y_ion_col, lwd = 2)
            
            text((yion_xlarge+PSManno_yion)/2,
                 rep(peptide_height - len_annoSpace, length(PSManno_yion)),
                 PSManno_ysmall, cex = 1, adj = c(0.5, -0.3),col = y_ion_col)
        }
    }
    
}
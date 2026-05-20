#' @rdname validatePSM
#'
#' @title Validate a PSM
#'
#' @description
#' Validate a PSM by checking some of its spectral features. Calling
#' `validatePSM()` applies the different checks mentioned below. Failing a check
#' does not guarantee a bad identification per se, but it does provide
#' information to assess the confidence in dubious PSMs or confirm already good
#' identifications.
#'
#' @param x A `Spectra` object with identifications. The identification
#' sequences need to be stocked in a variable name called with the parameter
#' `peptideVariable`. Refer to `Spectra::joinSpectraData` for more
#' informations.
#'
#' @param peptideVariable A `character` of length 1L, representing the peptide
#' sequence for that identification.
#'
#' @param fdr A `character` of length 1L, representing the FDR value for that
#' identification.
#'
#' @param ... Arguments passed down to `checkParentIonIntensity()`
#'
#' @returns A `data.frame()` that checks the different validation metrics.
#'
#' @importFrom Spectra spectraData peaksData
#'
#' @export
#'
#' @examples
#'
#' library(Spectra)
#'
#' data("psmBoekweg")
#' data("spBoekweg")
#'
#' ## Make sure you can join both using `joinSpectraData`:
#' head(psmBoekweg$pkey <- paste0(basename(psmBoekweg$filename),
#'                                sub("^.+scan=", "::", psmBoekweg$scannr)))
#' head(spBoekweg$pkey <- paste0(basename(spBoekweg$dataOrigin),
#'                               sub("^.+scan=", "::", spBoekweg$spectrumId)))
#'
#' sp <- Spectra::joinSpectraData(spBoekweg, psmBoekweg, by.x = "pkey")
#'
#' ## Add carbamidomethylation or other modifications if need be
#' ## See ?PTMods::addFixedModifications
#' (seq <- psmVariables(psmBoekweg)[["peptide"]])
#' (fdr <- psmVariables(psmBoekweg)[["fdr"]])
#' head(sp[["modSequences"]] <- addFixedModifications(sp[[seq]]))
#'
#' validatePSM(sp[1:20], peptideVariable = "modSequences", fdr = fdr)
#'
validatePSM <- function(x, peptideVariable = "peptide", fdr = "fdr", ...) {

    stopifnot(requireNamespace("Spectra"))
    stopifnot(inherits(x, "Spectra"))

    v <- Spectra::peaksData(x) ## all spectra (MS1 & MS2) peaks data

    sequenceIds <- Spectra::spectraData(x, peptideVariable)[, 1L]
    identifications <- which(!is.na(sequenceIds))

    x$sequence <- sequenceIds ## Create 'sequence' variable for labelFragments

    x_sub <- x[identifications] ## Fetch only PSMs
    sub_v <- v[identifications] ## Identified spectra peak data
    stripped_sequence <- PTMods::getCanonicalSequence(sequenceIds[identifications])

    stopifnot(length(sub_v) != 0)

    frags <- suppressWarnings(labelFragments(x_sub,
        type = c("a", "b", "c", "x", "y", "z")))
    frags_mz <- suppressWarnings(labelFragments(x_sub,
        type = c("b", "y"), what = "mz"))

    ab <- checkABpresence(x_sub, fragments = frags)

    xy <- checkXYpresence(x_sub, fragments = frags)

    overlap <- checkOverlap(x_sub, peptide = peptideVariable,
        fragments = frags, strippedSeq = stripped_sequence)

    shifts <- checkShiftConsistency(x_sub, fragments = frags)

    parent <- checkParentIonIntensity(x_sub, peaks = sub_v, ...)

    purity_ints <- checkPrecursorPurity(x)[identifications]

    res <- data.frame(spectrumId = Spectra::spectraData(x_sub, "spectrumId")[, 1L],
                      scanIndex = Spectra::scanIndex(x_sub),
                      peptide = sequenceIds[identifications],
                      canonicalSeq = stripped_sequence,
                      fdr = Spectra::spectraData(x_sub, fdr)[, 1L],
                      a2b2 = ab,
                      x2y2 = xy,
                      byOverlap = overlap,
                      shiftConsistency = shifts,
                      parentIonInt = parent,
                      precursorPurity = purity_ints)

    rownames(res) <- NULL
    return(res)
}

#' @rdname validatePSM
#'
#' @param fragments The result of `labelFragments()` on `x`.
#'
#' @returns `checkABpresence()` : `TRUE` if the a2-b2 fragments are both present.
#'
#' @examples
#'
#' library(Spectra)
#'
#' ## Build a minimal spectrum containing a- and b-ions for "PEPTIDE"
#' seq <- "PEPTIDE"
#' frags_ab <- calculateFragments(seq, type = c("a", "b"))
#' sp_ab <- DataFrame(msLevel = 2L, rtime = 100, sequence = seq)
#' sp_ab$mz <- list(sort(frags_ab$mz))
#' sp_ab$intensity <- list(rep(1e4, nrow(frags_ab)))
#' sp_ab <- Spectra(sp_ab)
#' ## Both a2 and b2 are matched: returns TRUE
#' checkABpresence(sp_ab)
#'
#' @export
checkABpresence <- function(x, fragments = NULL) {

    if (!length(fragments)) {
        labels <- labelFragments(x, type = c("a", "b"))
    } else {
        labels <- fragments
    }
    .abCouple <- function(labs) "a2" %in% labs & "b2" %in% labs
    unlist(lapply(labels, .abCouple))
}

#' @rdname validatePSM
#'
#' @param fragments The result of `labelFragments()` on `x`.
#'
#' @returns `checkXYpresence()` : `TRUE` if the x2-y2 fragments are both present.
#'
#' @examples
#'
#' ## Pass a fragment label list directly to checkXYpresence
#' sp_xy <- DataFrame(msLevel = 2L, rtime = 100)
#' sp_xy$mz <- list(c(100, 200, 300, 400))
#' sp_xy$intensity <- list(rep(1e4, 4))
#' sp_xy <- Spectra(sp_xy)
#' ## Fragment list containing x2 and y2: returns TRUE
#' checkXYpresence(sp_xy, fragments = list(c("x2", "y2", "y3")))
#' ## Fragment list missing x2: returns FALSE
#' checkXYpresence(sp_xy, fragments = list(c("y2", "y3", "y4")))
#'
#' @export
checkXYpresence <- function(x, fragments = NULL) {

    if (!length(fragments)) {
        labels <- labelFragments(x, type = c("x", "y"))
    } else {
        labels <- fragments
    }
    .xyCouple <- function(labs) "x2" %in% labs & "y2" %in% labs
    unlist(lapply(labels, .xyCouple))
}

#' @rdname validatePSM
#'
#' @param fragments The result of `labelFragments()` on `x`.
#'
#' @returns `checkOverlap()` : Detects a gap in the coverage of b- and y-ions.
#' Returns `FALSE` when a gap is found (b- and y-ions do not jointly cover the
#' full sequence), which may indicate an unsearched modification. Returns
#' `TRUE` when b- and y-ions overlap, further confirming the identification.
#'
#' @examples
#'
#' ## Sparse fragment coverage → gap detected (FALSE)
#' seq <- "PEPTIDE"
#' frags_all <- calculateFragments(seq, type = c("b", "y"))
#' frags_sparse <- frags_all[frags_all$ion %in% c("b1", "b2", "y1", "y2"), ]
#' sp_gap <- DataFrame(msLevel = 2L, rtime = 100, sequence = seq)
#' sp_gap$mz <- list(sort(frags_sparse$mz))
#' sp_gap$intensity <- list(rep(1e4, nrow(frags_sparse)))
#' sp_gap <- Spectra(sp_gap)
#' labs_gap <- suppressWarnings(labelFragments(sp_gap, type = c("b", "y")))
#' checkOverlap(sp_gap, fragments = labs_gap, strippedSeq = seq)
#'
#' ## Full b and y coverage → no gap (TRUE)
#' sp_full <- DataFrame(msLevel = 2L, rtime = 100, sequence = seq)
#' sp_full$mz <- list(sort(frags_all$mz))
#' sp_full$intensity <- list(rep(1e4, nrow(frags_all)))
#' sp_full <- Spectra(sp_full)
#' labs_full <- suppressWarnings(labelFragments(sp_full, type = c("b", "y")))
#' checkOverlap(sp_full, fragments = labs_full, strippedSeq = seq)
#'
#' @export
checkOverlap <- function(x, peptideVariable = "peptide",
                         fragments = NULL, strippedSeq = NULL) {

    if (length(strippedSeq)) {
        stripped_sequence <- strippedSeq
    } else {
        sequence_ids  <- Spectra::spectraData(x, peptideVariable)[, 1L]
        stripped_sequence <- PTMods::getCanonicalSequence(sequence_ids)
    }

    if (!length(fragments)) {
        x$sequence <- Spectra::spectraData(x, peptideVariable)[, 1L]
        labels <- suppressWarnings(labelFragments(x))
    } else {
        labels <- fragments
    }

    ## fetch all b and y fragments based on calculateFragments()
    ## including neutral losses
    index <- lapply(labels, function(x) as.integer(gsub("\\D", "", x)))
    b_ions <- lapply(labels, function(x) which(grepl("b", x)))
    y_ions <- lapply(labels, function(x) which(grepl("y", x)))

    ans <- vector(length = length(x))

    ## fetch only the highest fragments for both b and y fragments
    for (i in seq_along(x)) {
        max_b <- max(index[[i]][b_ions[[i]]], 0)
        max_y <- max(index[[i]][y_ions[[i]]], 0)
        ## n - max(y) > max(b), if respected = no overlap
        ans[i] <- nchar(stripped_sequence[i]) - max_y > max_b
    }
    return(!ans)
}

#' @rdname validatePSM
#'
#' @param fragments The result of `labelFragments()` on `x`.
#'
#' @returns `checkShiftConsistency()` : In case of modifications present:
#' returns the percentage of potential mass shifts actually matched. If not
#' applicable because there are no modifications: `NA`.
#'
#' @examples
#'
#' ## Unmodified sequence → no shift to assess, returns NA
#' seq <- "PEPTIDE"
#' frags_sc <- calculateFragments(seq)
#' sp_sc <- DataFrame(msLevel = 2L, rtime = 100, sequence = seq)
#' sp_sc$mz <- list(sort(frags_sc$mz))
#' sp_sc$intensity <- list(rep(1e4, nrow(frags_sc)))
#' sp_sc <- Spectra(sp_sc)
#' labs_sc <- suppressWarnings(labelFragments(sp_sc))
#' checkShiftConsistency(sp_sc, fragments = labs_sc)
#'
#' @export
checkShiftConsistency <- function(x, fragments = NULL, strippedSeq = NULL) {

    if (length(fragments)) {
        labels = fragments
    } else {
        x$sequence <- Spectra::spectraData(x, peptideVariable)[, 1L]
        labels <- suppressWarnings(labelFragments(x))
    }

    if (length(strippedSeq)) {
        stripped_sequence <- strippedSeq
    } else {
        stripped_sequence <- PTMods::getCanonicalSequence(x$sequence)
    }

    index <- lapply(labels, function(x) as.integer(gsub("\\D", "", x)))
    b_ions <- lapply(labels, function(x) which(grepl("b", x)))
    y_ions <- lapply(labels, function(x) which(grepl("y", x)))

    ans <- vector(length = length(x))

    for (i in seq_along(x)) {

        if (is.na(x$sequence[i])) {
            ans[i] <- NA
        } else if (x$sequence[i] == stripped_sequence[i]) {
            ans[i] <- NA
        } else {
            parsed_mods <- PTMods:::.parseModifiedSequence(x$sequence[i])
            pep_len <- nchar(stripped_sequence[i])
            mod_pos <- which(parsed_mods != 0)
            max_b <- mod_pos[1] - 1
            max_y <- pep_len - mod_pos[length(mod_pos)]

            b_matches <- sum(unique(index[[i]][b_ions[[i]]]) > max_b)
            y_matches <- sum(unique(index[[i]][y_ions[[i]]]) > max_y)

            if (sum(max_b, max_y) != 0) {
                ans[i] <- sum(b_matches, y_matches)/sum(max_b, max_y)
                } else ans[i] <- 1
        }
    }
    return(ans)
}

#' @rdname validatePSM
#'
#' @param peaks The spectrum peak data (result from `Spectra::peaksData()`).
#'
#' @param tolerance `Numeric(1L)` The tolerance to use when matching peaks.
#'
#' @param ppm `Numeric(1L)` The ppm value to use when matching peaks that is
#'   added to `tolerance`.
#'
#' @importFrom MsCoreUtils closest
#'
#' @returns `checkParentIonIntensity()` : The relative intensity of the parent
#' ion over the base peak. A value close to 1 indicates poor fragmentation.
#'
#' @examples
#'
#' ## Precursor ion is the most intense peak → ratio of 1
#' sp_p <- DataFrame(msLevel = 2L, rtime = 100, precursorMz = 500.0)
#' sp_p$mz <- list(c(100, 200, 300, 400, 500))
#' sp_p$intensity <- list(c(100, 100, 100, 100, 1e6))
#' sp_p <- Spectra(sp_p)
#' checkParentIonIntensity(sp_p)
#'
#' ## Precursor absent from spectrum → ratio of 0
#' sp_np <- DataFrame(msLevel = 2L, rtime = 100, precursorMz = 999.0)
#' sp_np$mz <- list(c(100, 200, 300, 400, 500))
#' sp_np$intensity <- list(c(100, 100, 100, 100, 1e6))
#' sp_np <- Spectra(sp_np)
#' checkParentIonIntensity(sp_np)
#'
#' @export
checkParentIonIntensity <- function (x,
                                     peaks = NULL,
                                     tolerance = 0,
                                     ppm = 20) {

    stopifnot(requireNamespace("Spectra"))
    stopifnot(inherits(x, "Spectra"))

    k <- numeric()
    if (length(peaks)) {
        v <- peaks
    } else {
        v <- Spectra::peaksData(x)
    }
    precursors <- Spectra::precursorMz(x)

    for (i in seq_along(x)) {

        peak <- v[[i]]
        precursor_index <- MsCoreUtils::closest(precursors[i],
                                                peak[, "mz"],
                                                tolerance = tolerance,
                                                ppm = ppm)
        if (!is.na(precursor_index)) {
            ints <- peak[precursor_index, "intensity"]
        } else {ints <- 0}
        max_ints <- max(peak[, "intensity"])
        k[i] <- ints/max_ints
    }
    return(k)
}

#' @rdname validatePSM
#'
#' @examples
#'
#' ## Check precursor purity on the bundled spectra dataset
#' data("spBoekweg")
#' purity <- checkPrecursorPurity(spBoekweg[1:5])
#' purity
#'
#' @param x A `Spectra` object containing **both MS1 and MS2 spectra** from
#'   the same run(s). MS1 spectra are used as the source for isolation window
#'   peak data; each MS2 is paired with the nearest preceding MS1 scan (sorted
#'   by retention time within each `dataOrigin`).
#'
#' @param tolerance `Numeric(1)` Absolute m/z half-width (in Da) of the
#'   isolation window used when `useReportedIsolationWindow = FALSE` (default
#'   `0.05`).
#'
#' @param ppm `Numeric(1)` Additional m/z-proportional tolerance (in ppm)
#'   added to `tolerance` when defining the isolation window
#'   (`default 0`).
#'
#' @param useReportedIsolationWindow `logical(1)` If `TRUE`, use the
#'   `isolationWindowLowerMz` / `isolationWindowUpperMz` metadata stored in the
#'   spectra rather than computing the window from `tolerance` and `ppm`.
#'   Defaults to `FALSE`.
#'
#' @param BPPARAM A `BiocParallelParam` instance controlling parallel
#'   evaluation (one job per unique `dataOrigin`). Defaults to
#'   `BiocParallel::SerialParam()`.
#'
#' @importFrom Spectra precursorPurity
#'
#' @returns `checkPrecursorPurity()` returns a `numeric` vector of length
#'   `length(x)`. Each value is the ratio of the most-intense peak to the
#'   total intensity within the isolation window of the corresponding MS1 scan,
#'   as computed by `Spectra::precursorPurity()`. MS1 spectra and MS2 spectra
#'   with no preceding MS1 scan return `NA`.
#'
#' @export
checkPrecursorPurity <- function(x, tolerance = 0.05, ppm = 0,
                                 useReportedIsolationWindow = FALSE,
                                 BPPARAM = BiocParallel::SerialParam()) {

    stopifnot(inherits(x, "Spectra"))

    Spectra::precursorPurity(x,
                            tolerance = tolerance,
                            ppm = ppm,
                            useReportedIsolationWindow =
                                useReportedIsolationWindow,
                            BPPARAM = BPPARAM)
}
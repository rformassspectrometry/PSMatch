##' Example peptide-spectrum match data
##'
##' @description
##'
##' `psmBoekweg` is a `PSM` object containing 10,000 peptide-spectrum
##' matches randomly sampled from a larger dataset. The PSMs were
##' identified using the Sage search engine.
##'
##' @format
##'
##' A `PSM` object with 10,000 rows. Each PSM contains the following
##' metadata columns (among others):
##'
##' \describe{
##'   \item{`psm_id`}{Unique identifier for the PSM}
##'   \item{`peptide`}{Peptide sequence}
##'   \item{`proteins`}{Protein identifiers}
##'   \item{`filename`}{Source mzML file name}
##'   \item{`scannr`}{Scan number}
##'   \item{`rank`}{Rank of the PSM}
##'   \item{`charge`}{Precursor charge state}
##'   \item{`hyperscore`}{Sage hyperscore}
##'   \item{`peptide_q`}{Peptide-level q-value (FDR)}
##'   \item{`protein_q`}{Protein-level q-value (FDR)}
##'   \item{`rt`}{Retention time (seconds)}
##' }
##'
##' @details
##'
##' This is an example subset of the data available through
##' `?MsDataHub::Boekweg2022()`.
##'
##' @seealso
##'
##' * [spBoekweg] for the corresponding `Spectra` object.
##' * [PSM()] for the PSM class constructor.
##' * The `inst/scripts/make_BoekwegData.R` script documenting data
##'   creation.
##'
##' @examples
##' data(psmBoekweg)
##' psmBoekweg
##'
##' ## Number of PSMs
##' length(psmBoekweg)
##'
##' ## Access peptide sequences
##' head(psmBoekweg$peptide)
##'
##' ## Filter to rank 1 PSMs only
##' psmRank1 <- filterPsmRank(psmBoekweg)
##' length(psmRank1)
##'
##' ## Summary of q-values
##' summary(psmBoekweg$peptide_q)
"psmBoekweg"


##' Example mass spectrometry data
##'
##' @description
##'
##' `spBoekweg` is a `Spectra` object containing 4,766 mass spectra
##' (MS1 and MS2) that are linked to peptide identifications in the
##' `psmBoekweg` dataset.
##'
##' @format
##'
##' A `Spectra` object with 4,766 spectra stored in memory using
##' `MsBackendMemory()`. The spectra contain standard metadata
##' variables including:
##'
##' \describe{
##'   \item{`msLevel`}{MS level (1 or 2)}
##'   \item{`rtime`}{Retention time (seconds)}
##'   \item{`precursorMz`}{Precursor m/z (MS2 only)}
##'   \item{`precursorCharge`}{Precursor charge (MS2 only)}
##'   \item{`dataOrigin`}{Source file name}
##'   \item{`spectrumId`}{Spectrum identifier}
##' }
##'
##' @details
##'
##' This is an example subset of the data available through
##' `?MsDataHub::Boekweg2022()`.
##'
##' The spectra can be linked to the `psmBoekweg` dataset using a
##' common key based on the file name and scan number. See the vignette
##' `validatePSM` for an example.
##'
##' @seealso
##'
##' * [psmBoekweg] for the corresponding `PSM` object with peptide
##'   identifications.
##' * [Spectra::Spectra()] for the Spectra class from the Spectra
##'   package.
##' * The `inst/scripts/make_BoekwegData.R` script documenting data
##'   creation.
##'
##' @examples
##' data(spBoekweg)
##' spBoekweg
##'
##' ## Number of spectra
##' length(spBoekweg)
##'
##' ## MS levels
##' table(spBoekweg$msLevel)
##'
##' ## Retention time range (in seconds)
##' range(spBoekweg$rtime)
##'
##' ## Plot base peak chromatogram
##' plot(spBoekweg$rtime, spBoekweg$basePeakIntensity,
##'      xlab = "Retention time (s)", ylab = "Base peak intensity",
##'      type = "l")
##'
##' ## Link to PSM data using a common key
##' data(psmBoekweg)
##' psmBoekweg$pkey <- paste0(
##'     basename(psmBoekweg$filename),
##'     sub("^.+scan=", "::", psmBoekweg$scannr)
##' )
##' spBoekweg$pkey <- paste0(
##'     basename(spBoekweg$dataOrigin),
##'     sub("^.+scan=", "::", spBoekweg$spectrumId)
##' )
##'
##' ## Find spectra with identifications
##' sum(spBoekweg$pkey %in% psmBoekweg$pkey)
"spBoekweg"

#' PSMatch: Handling and Managing Peptide Spectrum Matches
#'
#' The PSMatch package offers functionality to load, manage and
#' analyse Peptide Spectrum Matches as generated in mass
#' spectrometry-based proteomics. The three main objects and concepts
#' that are proposed in this package are described below, and are
#' aimed to proteomics practitioners to explore and understand their
#' identification data better.
#'
#' @section PSM objects:
#'
#' As mentioned in the [PSM()] manual page, The `PSM` class is a
#' simple class to store and manipulate peptide-spectrum matches. The
#' class encapsulates PSM data as a DataFrame (or more specifically a
#' `DFrame`) with additional lightweight metadata annotation. PSM
#' objects are typically creatd from XML-based mzID files or
#' `data.frames` imported from spreadsheets. It is them possible to
#' apply widely used filters (such as removal of decoy hits, PSMs of
#' rank > 1, ...) as described in [filterPSMs()].
#'
#' @section Adjacency matrices:
#'
#' PSM data, as produced by all proteomics search engines, is exported
#' as a table like structure where PSM are documented along the rows
#' by variables such as identification scores, peptides sequences,
#' modifications and the protein which that peptide comes from. There
#' is always a level of ambiguity in such data, as peptides can be
#' mapped to mutliple proteins; they are then called shared peptides,
#' as opposed to unique peptides.
#'
#' One convenient way to store the relation between peptides and
#' proteins is as a peptide-by-protein adjacency matrix. Such matrices
#' can be generated from PSM object or vectors using the
#' [makeAdjacencyMatrix()].
#'
#' The [describePeptides()] and [describeProteins()] functions are
#' also helpful is tally the number of unique and shared peptides and
#' the number of proteins composed of unique or shared peptides, or a
#' combination thereof.
#'
#' @section Connected Components:
#'
#' Once we model the peptide-to-protein relations explicitly, it
#' becomes possible to perform computations on the proteins that are
#' grouped by the peptides they share. These groups are mathematically
#' defined as connected components, which are implemented as
#' [ConnectedComponents()] objects.
#'
#' @section Fragment ions:
#'
#' The package also provides functionality to calculate ions produced
#' by the fragmentation of a peptides (see [calculateFragments()]) and
#' annotated MS2 [Spectra::Spectra()] objects (see [addFragments()]).
#'
#' @docType package
#'
#' @name PSMatch
NULL

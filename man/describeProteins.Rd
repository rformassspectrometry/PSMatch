% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/describe.R
\name{describeProteins}
\alias{describeProteins}
\alias{describePeptides}
\title{Describe protein and peptide compositions}
\usage{
describeProteins(object, ...)

describePeptides(object, ...)
}
\arguments{
\item{object}{Either an instance of class \code{Matrix}, \code{\link[=PSM]{PSM()}} or
\code{\link[=ConnectedComponents]{ConnectedComponents()}}.}

\item{...}{Additional arguments passed to \code{\link[=makeAdjacencyMatrix]{makeAdjacencyMatrix()}}.}
}
\value{
\code{describePeptides()} invisibly return the table of unique
and shared peptides. \code{describeProteins()} invisibly returns a
\code{data.frame} with logicals indicating the unique/shared
peptide composition of proteins. Both functions are used for
their side effects of printing a short descriptive output
about peptides and proteins.
}
\description{
It is important to explore PSM results prior to any further
downstream analysies. Two functions, that work on \code{\link[=PSM]{PSM()}} and
\code{\link[=ConnectedComponents]{ConnectedComponents()}} objects can be used for this:
\itemize{
\item The \code{describeProteins()} function describe protein composition
in terms of unique and shared peptides.
\item The \code{describePeptides()} function describe unique/shared peptide
composition.
}
}
\examples{
f <- msdata::ident(full.names = TRUE, pattern = "TMT")
basename(f)
psm <- PSM(f) |>
       filterPsmDecoy() |>
       filterPsmRank()

describePeptides(psm)
describeProteins(psm)
}

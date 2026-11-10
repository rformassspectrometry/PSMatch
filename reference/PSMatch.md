# PSMatch: Handling and Managing Peptide Spectrum Matches

The PSMatch package offers functionality to load, manage and analyse
Peptide Spectrum Matches as generated in mass spectrometry-based
proteomics. The four main objects and concepts that are proposed in this
package are described below, and are aimed to proteomics practitioners
to explore and understand their identification data better.

## PSM objects

As mentioned in the
[`PSM()`](https://rformassspectrometry.github.io/PSMatch/reference/PSM.md)
manual page, The `PSM` class is a simple class to store and manipulate
peptide-spectrum matches. The class encapsulates PSM data as a DataFrame
(or more specifically a `DFrame`) with additional lightweight metadata
annotation. PSM objects are typically creatd from XML-based mzID files
or `data.frames` imported from spreadsheets. It is then possible to
apply widely used filters (such as removal of decoy hits, PSMs of rank
\> 1, ...) as described in
[`filterPSMs()`](https://rformassspectrometry.github.io/PSMatch/reference/filterPSMs.md).

## Adjacency matrices

PSM data, as produced by all proteomics search engines, is exported as a
table-like structure where PSM are documented along the rows by
variables such as identification scores, peptides sequences,
modifications and the protein which the peptides originate from. There
is always a level of ambiguity in such data, as peptides can be mapped
to mutliple proteins; they are then called shared peptides, as opposed
to unique peptides.

One convenient way to store the relation between peptides and proteins
is as a peptide-by-protein adjacency matrix. Such matrices can be
generated from PSM object or vectors using the
[`makeAdjacencyMatrix()`](https://rformassspectrometry.github.io/PSMatch/reference/adjacencyMatrix.md)
function.

The
[`describePeptides()`](https://rformassspectrometry.github.io/PSMatch/reference/describeProteins.md)
and
[`describeProteins()`](https://rformassspectrometry.github.io/PSMatch/reference/describeProteins.md)
functions are also helpful to tally the number of unique and shared
peptides and the number of proteins composed of unique or shared
peptides, or a combination thereof.

## Connected Components

Once we model the peptide-to-protein relations explicitly using an
adjacency matrix, it becomes possible to perform computations on the
proteins that are grouped by the peptides they share. These groups are
mathematically defined as connected components, which are implemented as
[`ConnectedComponents()`](https://rformassspectrometry.github.io/PSMatch/reference/ConnectedComponents.md)
objects.

## Fragment ions

The package also provides functionality to calculate ions produced by
the fragmentation of a peptides (see
[`calculateFragments()`](https://rformassspectrometry.github.io/PSMatch/reference/calculateFragments.md))
and annotated MS2
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
objects (see
[`labelFragments()`](https://rformassspectrometry.github.io/PSMatch/reference/labelFragments.md)).

## Vignettes

A couple of vignette describe how to several of these concepts through
illustrative use-cases. Use `vignette(package = "PSMatch")` to get a
list and open them directly in `R` or read them online on the package's
[webpage](https://rformassspectrometry.github.io/PSMatch/).

## See also

Useful links:

- <https://github.com/RforMassSpectrometry/PSM>

- Report bugs at <https://github.com/RforMassSpectrometry/PSM/issues>

## Author

**Maintainer**: Laurent Gatto <laurent.gatto@uclouvain.be>
([ORCID](https://orcid.org/0000-0002-1520-2268))

Authors:

- Laurent Gatto <laurent.gatto@uclouvain.be>
  ([ORCID](https://orcid.org/0000-0002-1520-2268))

- Johannes Rainer <Johannes.Rainer@eurac.edu>
  ([ORCID](https://orcid.org/0000-0002-6977-7147))

- Sebastian Gibb <mail@sebastiangibb.de>
  ([ORCID](https://orcid.org/0000-0001-7406-4443))

Other contributors:

- Samuel Wieczorek <samuel.wieczorek@cea.fr> \[contributor\]

- Thomas Burger <thomas.burger@cea.fr> \[contributor\]

- Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
  ([ORCID](https://orcid.org/0009-0008-1257-2416)) \[contributor\]

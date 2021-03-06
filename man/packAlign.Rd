% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/packAlign.R
\name{packAlign}
\alias{packAlign}
\title{Global Alignment with VSEARCH}
\usage{
packAlign(
  packMatches,
  Genome,
  identity = 0,
  threads = 1,
  identityDefinition = 2,
  maxWildcards = 0.05,
  saveFolder,
  vSearchPath = "vsearch"
)
}
\arguments{
\item{packMatches}{A dataframe of potential Pack-TYPE transposable elements, 
in the format given by \code{\link{packSearch}}. This 
dataframe is in the format produced by coercing a 
\code{link[GenomicRanges:GRanges-class]{GRanges}} 
object to a dataframe: \code{data.frame(GRanges)}. 
Will be saved as a FASTA file for VSEARCH.}

\item{Genome}{A DNAStringSet object containing sequences referred to in 
\code{packMatches} (the object originally used to predict 
the transposons \code{\link{packSearch}}).}

\item{identity}{The sequence identity of two transposable elements in 
\code{packMatches} required to be grouped into a cluster.}

\item{threads}{The number of threads to be used by VSEARCH.}

\item{identityDefinition}{The pairwise identity definition used by VSEARCH. 
Defaults to 2, the standard VSEARCH definition.}

\item{maxWildcards}{The maximal allowable proportion of wildcards in the 
sequence of each match (defaults to \code{0.05}).}

\item{saveFolder}{The folder to save saveFolder files (uc, blast6out, FASTA)}

\item{vSearchPath}{When the package is run on windows systems, the 
location of the VSEARCH executable file must be 
given; this should be left as default on 
Linux/MacOS systems.}
}
\value{
Saves alignment information, including a \code{uc}, 
\code{blast6out} and a pairwise alignment \code{fasta} 
file, to the specified location. Returns the uc summary 
file generated by the alignment.
}
\description{
A global pairwise alignment of pack-TYPE elements by 
sequence similarity. mIt may be useful to run 
\code{\link{packClust}} to identify groups of similar 
transposable elements, before generating alignments of 
each group.
}
\note{
In order to align sequences using VSEARCH, the 
executable file must first be installed.
}
\examples{
data(arabidopsisThalianaRefseq)
data(packMatches)

# packAlign run on a Linux/MacOS system
\dontrun{
    packAlign(packMatches, Genome)
}

# packAlign run on a Windows system
\dontrun{
    packAlign(packMatches, Genome, 
            vSearchPath = "path/to/vsearch/vsearch.exe")
}

}
\references{
VSEARCH may be downloaded from 
\url{https://github.com/torognes/vsearch}, along with a 
manual documenting the program's parameters. See 
\url{https://www.ncbi.nlm.nih.gov/pubmed/27781170} for 
further information.
}
\seealso{
\code{\link{tirClust}}, \code{\link{packClust}}, 
\code{\link{readBlast}}, \code{\link{readUc}},
\code{\link{filterWildcards}}, \code{\link{packSearch}}
}
\author{
Jack Gisby
}

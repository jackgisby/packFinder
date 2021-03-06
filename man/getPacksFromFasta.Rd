% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPacksFromFasta.R
\name{getPacksFromFasta}
\alias{getPacksFromFasta}
\title{Retrieve Saved packFinder Results (.fasta)}
\usage{
getPacksFromFasta(file)
}
\arguments{
\item{file}{Path to predicted transposons in FASTA format.}
}
\value{
Dataframe in the format used by \code{\link{packSearch}}.
}
\description{
Retrieves a dataframe of potential Pack-TYPE elements, 
previously saved using \code{\link{packSearch}} followed 
by \code{\link{packsToFasta}}. Parses the .fasta file 
and title field containing:
\itemize{
    \item seqnames - name of origin sequence
    \item start - transposon base start position on 
    origin sequence
    \item end - transposon base end position on origin 
    sequence
    \item width - width of transposon
    \item strand - direction of transposon 
    ("+", "-" or "*")
    \item TSD - terminal site duplication (TSD) sequence
}
}
\examples{
data(arabidopsisThalianaRefseq)
data(packMatches)

packMatches <- getPacksFromFasta(
    system.file("extdata", "packMatches.fasta", package = "packFinder")
)

}
\seealso{
\code{\link{packsToFasta}}, \code{\link{packSearch}}
}
\author{
Jack Gisby
}

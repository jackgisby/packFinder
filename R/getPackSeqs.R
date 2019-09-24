#' @title
#' Extract Sequences of Pack-TYPE Elements
#'
#' @description
#' Method to quickly extract the sequences of predicted Pack-TYPE elements (as
#' created by \code{\link{packSearch}}).
#'
#' @param packMatches
#' A dataframe containing genomic ranges and names referring
#' to sequences to be extracted. Can be obtained from \code{\link{packSearch}}
#' or generated from a \code{\link[GenomicRanges]{GRanges}} object, after
#' conversion to a dataframe. Must contain the following features:
#' \itemize {
#'   \item start - the predicted element's start base sequence position.
#'   \item end - the predicted element's end base sequence position.
#'   \item seqnames - character string referring to the sequence name in
#'   \code{Genome} to which \code{start} and \code{end} refer to.
#' }
#'
#' @param Genome
#' A DNAStringSet object containing sequences referred to in \code{packMatches}
#' (the object originally used to predict the transposons
#' \code{\link{packSearch}}).
#'
#' @param output
#' The type of object to be returned:
#' \itemize {
#'  \item output = "DNAStringSet", returns a
#'  \code{\link[Biostrings]{DNAStringSet}} object (default).
#'  \item output = "character", returns a \code{character} vector.
#'  }
#'
#' @examples \dontrun{
#'
#' data(arabidopsisThalianaRefseqSubset)
#'
#' packMatches <- packSearch(
#'   Biostrings::DNAString("CACTACAA"),
#'   arabidopsisThalianaRefseq,
#'   elementLength = c(300, 3500),
#'   tsdLength = 3)
#'
#' packSeqs <- getPackSeqs(packMatches,
#'                         arabidopsisThalianaRefseq)}
#'
#' @author
#' Jack Gisby
#'
#' @return
#' The transposon sequences extracted from \code{packMatches}. At default
#' returns the sequences as a \code{\link[Biostrings]{DNAStringSet}} or, if
#' \code{output} is set to "character", returns a character vector.
#'
#' @export

getPackSeqs <- function(packMatches,
                        Genome,
                        output = "DNAStringSet") {
  if (output != "string" & output != "DNAStringSet") {
    stop("Argument 'output' must be specified as 'string' or 'DNAStringSet'")
  }

  seqs <- mapply(function(start,
                            end,
                            seqnames,
                            Genome) {
    return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][start:end]))
  },
  packMatches$start,
  packMatches$end,
  packMatches$seqnames,
  MoreArgs = list(Genome = Genome)
  )

  if (output == "string") {
    return(seqs)
  } else if (output == "DNAStringSet") {
    seqs <- Biostrings::DNAStringSet(seqs)
    seqs@ranges@NAMES <- as.character(packMatches$seqnames)
    return(seqs)
  }
}

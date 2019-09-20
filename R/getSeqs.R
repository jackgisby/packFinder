#' @title Get Sequences from Ranges.
#' @description Gets the sequences referred to in \code{packMatches} and returns the sequences added as an additional column to the dataframe..
#' @param packMatches A dataframe containing genomic ranges and names referring to sequences to be extracted.
#' @param Genome A DNAStringSet object containing sequences referred to in \code{packMatches}
#' @author Jack Gisby
#' @details
#' Function intended for internal use.
#' @return The dataframe \code{packMatches} with an additional \code{seq} feature containing extracted sequences as characters.

getSeqs <- function(packMatches,
                    Genome) {

  packMatches %>%
    dplyr::mutate(seq = mapply(function(start, end, seqnames, Genome) {
      return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][start:end]))},
      start,
      end,
      seqnames,
      MoreArgs = list(Genome = Genome))) %>%
    return()
}

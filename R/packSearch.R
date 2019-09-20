#' @title packFind Algorithm Pipeline
#' @description General use pipeline function for the Pack-TYPE transposon finding algorithm.
#' @param subSeq A DNAString object containing the TIR sequence to be searched for.
#' @param Genome A DNAStringSet object to be searched.
#' @param mismatch The maximum edit distance to be considered for TIR matches (indels + substitions). See \code{\link[Biostrings]{matchPattern}}.
#' @param elementLength The maximum element length to be considered, as a vector of two integers. E.g. \code{c(300, 3500)}
#' @param tsdLength Integer referring to the length of the flanking TSD region.
#' @author Jack Gisby
#' @details
#' Finds potential pack-TYPE elements based on:
#' \itemize{
#'   \item Similarity of TIR sequence to \code{subSeq}
#'   \item Proximity of potential TIR sequences
#'   \item Directionality of TIR sequences
#'   \item Similarity of TSD sequences
#' }
#' @return A dataframe, \code{packMatches}, containing elements identified by the algorithm. These may be autonomous or pack-TYPE elements.
#' @export

packSearch <- function(subSeq,
                       Genome,
                       mismatch = 0,
                       elementLength,
                       tsdLength) {

  # perform initial search for TIR matches and get related TSD sequences
  print("Getting forward matches")
  forwardMatches <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), strand = character())
  forwardMatches <- identifyTirMatches(forwardMatches, subSeq, Genome, mismatch = mismatch, strand = "+")
  forwardMatches <- getTsds(forwardMatches, Genome, tsdLength, direction = "+")

  print("Getting reverse matches")
  reverseMatches <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), strand = character())
  reverseMatches <- identifyTirMatches(reverseMatches, Biostrings::reverseComplement(subSeq), Genome, mismatch = mismatch, strand = "-")
  reverseMatches <- getTsds(reverseMatches, Genome, tsdLength, direction = "-")

  # case: no matches
  if (length(forwardMatches[, 1]) == 0 | length(reverseMatches[, 1]) == 0) {
    print("No matches identified")
    return(NULL)
  }

  # determine potential transposable elements based on nearby elements and TSD sequences
  print("Filtering matches based on TSD sequences")
  packMatches <- identifyPotentialPackElements(forwardMatches, reverseMatches, Genome, elementLength)
  packMatches <- getTsds(packMatches, Genome, tsdLength, "+")
  packMatches <- getSeqs(packMatches, Genome)

  print("Initial filtering complete")
  return(packMatches)
}

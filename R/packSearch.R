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
#'
#' The algorithm finds potential forward and reverse TIR sequences using \code{\link{identifyTirMatches}} and their associated TSD sequence via \code{\link{getTsds}}. The main filtering stage, \code{\link{identifyPotentialPackElements}}, filters matches to obtain a dataframe of potential PACK elements. Note that this pipeline does not consider the possibility of discovered elements being autonomous elements, so it is recommended to cluster and/or BLAST elements for further analysis.
#' @return A dataframe, \code{packMatches}, containing elements identified by the algorithm. These may be autonomous or pack-TYPE elements.
#' @export

packSearch <- function(subSeq,
                       Genome,
                       mismatch = 0,
                       elementLength,
                       tsdLength) {

  # perform initial search for TIR matches and get related TSD sequences
  message("Getting forward matches")
  forwardMatches <- identifyTirMatches(
    subSeq = subSeq,
    Genome = Genome,
    mismatch = mismatch,
    strand = "+"
  )

  forwardMatches <- getTsds(
    tirMatches = forwardMatches,
    Genome = Genome,
    tsdLength = tsdLength,
    strand = "+"
  )

  message("Getting reverse matches")
  reverseMatches <- identifyTirMatches(
    subSeq = Biostrings::reverseComplement(subSeq),
    Genome = Genome,
    mismatch = mismatch,
    strand = "-"
  )

  reverseMatches <- getTsds(
    tirMatches = reverseMatches,
    Genome = Genome,
    tsdLength = tsdLength,
    strand = "-"
  )

  # case: no matches
  if (length(forwardMatches[, 1]) == 0 | length(reverseMatches[, 1]) == 0) {
    message("No matches identified")
    return(NULL)
  }

  # determine potential transposable elements based on nearby elements and TSD sequences
  message("Filtering matches based on TSD sequences")
  packMatches <- identifyPotentialPackElements(
    forwardMatches = forwardMatches,
    reverseMatches = reverseMatches,
    Genome = Genome,
    elementLength = elementLength
  )

  packMatches <- getTsds(
    tirMatches = packMatches,
    Genome = Genome,
    tsdLength = tsdLength,
    strand = "+"
  )

  message("Initial filtering complete")
  return(packMatches)
}

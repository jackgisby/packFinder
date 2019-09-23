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
  if (!is.integer(mismatch) | !is.integer(tsdLength)) {
    stop("Arguments 'mismatch' and 'tsdLength' must be integers")
  }

  if (!is.vector(elementLength) | length(elementLength) != 2) {
    stop("Argument 'elementLength' must be a vector of minimum and maximum transposon lengths")
  }

  if (!is.integer(elementLength[1]) | !is.integer(elementLength[2])) {
    stop("Vector 'elementLength' must contain integers")
  }

  if (typeof(subSeq) != "DNAString") {
    if (!is.character(subSeq)) {
      stop("Argument 'subSeq' must be of type Biostrings::DNAString or character")
    } else {
      subSeq <- Biostrings::DNAString(subSeq)
    }
  }

  if (typeof(Genome) != "DNAStringSet") {
    stop("Argument 'Genome' must be of type Biostrings::DNAStringSet.
         You may convert files using Biostrings::readDNAStringSet
         or convert objects using Biostrings::DNAStringSet")
  }

  # perform initial search for TIR matches and get related TSD sequences
  message("Getting forward matches")
  forwardMatches <- identifyTirMatches(
    subSeq = subSeq,
    Genome = Genome,
    mismatch = mismatch,
    strand = "+"
  )

  forwardMatches$TSD <- getTsds(
    tirMatches = forwardMatches,
    Genome = Genome,
    tsdLength = tsdLength,
    strand = "+"
  )

  message(length(forwardMatches[, 1]), " forward matches identified.")

  message("Getting reverse matches")
  reverseMatches <- identifyTirMatches(
    subSeq = Biostrings::reverseComplement(subSeq),
    Genome = Genome,
    mismatch = mismatch,
    strand = "-"
  )

  reverseMatches$TSD <- getTsds(
    tirMatches = reverseMatches,
    Genome = Genome,
    tsdLength = tsdLength,
    strand = "-"
  )

  message(length(reverseMatches[, 1]), " reverse matches identified.")

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

  packMatches$TSD <- getTsds(
    tirMatches = packMatches,
    Genome = Genome,
    tsdLength = tsdLength,
    strand = "+"
  )

  message("Initial filtering complete. ", length(packMatches[, 1]), " elements predicted.")
  return(packMatches)
}

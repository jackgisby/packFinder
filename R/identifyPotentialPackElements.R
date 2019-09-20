#' Primary filtering stage for the \code{packSearch} algorithm. Identifies potential packTYPE transposable elements based on proximity of matching inverted repeats and equality of TSD sequences.
#' @param forwardMatches A dataframe containing genomic ranges and names referring to forwards-facing TIR sequences and their respective TSD sequences.
#' @param reverseMatches A dataframe containing genomic ranges and names referring to forwards-facing TIR sequences and their respective TSD sequences.
#' @param Genome A DNAStringSet object containing the matches referred to in \code{forwardMatches} and \code{reverseMatches}
#' @param element.length A vector of two integers containing the minimum and maximum transposable element length.
#' @return A dataframe, \code{packMatches}, containing the locations of potential packTYPE transposable elements in \code{Genome}.
#' @export



identifyPotentialPackElements <- function(forwardMatches,
                                          reverseMatches,
                                          Genome,
                                          elementLength) {

  packMatches <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), strand = character())

  for(forwardMatch in 1:length(forwardMatches[,1])) {

    forwardRepeat <- forwardMatches[forwardMatch,]
    chr <- as.character(forwardRepeat[[1]])
    searchRange <- c(forwardRepeat$start + elementLength[1], forwardRepeat$start + elementLength[2])

    if(searchRange[2] > length(Genome[Genome@ranges@NAMES == chr][[1]])) {
      searchRange[2] <- length(Genome[Genome@ranges@NAMES == chr][[1]])
    }

    reverseRepeats <- filter(reverseMatches,
                             seqnames == as.character(forwardRepeat$seqnames) &
                               end > searchRange[1] &
                               end < searchRange[2] &
                               strand == "-"  &
                               TSD == as.character(forwardRepeat$TSD))

    if(length(reverseRepeats[,1]) > 0) {
      for(reverseMatch in 1:length(reverseRepeats[,1])) {
        packMatches <- rbind(packMatches, data.frame(seqnames = forwardRepeat$seqnames,
                                                     start = forwardRepeat$start,
                                                     end = reverseRepeats[reverseMatch,]$end,
                                                     width = reverseRepeats[reverseMatch,]$end - forwardRepeat$start,
                                                     strand = "*"))
      }
    }
  }
  return(packMatches)
}

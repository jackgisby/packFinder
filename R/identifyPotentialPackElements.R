#' @title Pack Element Filtering
#' @description Primary filtering stage for the \code{packSearch} algorithm. Identifies potential packTYPE transposable elements based on proximity of matching inverted repeats and equality of TSD sequences.
#' @param forwardMatches A dataframe containing genomic ranges and names referring to forwards-facing TIR sequences and their respective TSD sequences.
#' @param reverseMatches A dataframe containing genomic ranges and names referring to forwards-facing TIR sequences and their respective TSD sequences.
#' @param Genome A DNAStringSet object containing the matches referred to in \code{forwardMatches} and \code{reverseMatches}
#' @param elementLength A vector of two integers containing the minimum and maximum transposable element length.
#' @author Jack Gisby
#' @details
#' Used by \code{\link{packSearch}} as a primariy filtering stage. Identifies matches likely to be transposons based on their TIR region, from \code{\link{identifyTirMatches}}, and their TSD region, from \code{\link{getTsds}}. It is recommended to use the general pipeline function \code{\link{packSearch}} for identification of potential pack elements, however each stage may be called individually.
#' @return A dataframe, \code{packMatches}, containing the locations of potential packTYPE transposable elements in \code{Genome}.
#' @export



identifyPotentialPackElements <- function(forwardMatches,
                                          reverseMatches,
                                          Genome,
                                          elementLength) {
  packMatches <- data.frame(
    seqnames = character(),
    start = integer(),
    end = integer(),
    width = integer(),
    strand = factor()
  )

  # for each forward match
  for (forwardMatch in 1:length(forwardMatches[, 1])) {
    forwardRepeat <- forwardMatches[forwardMatch, ]
    chr <- as.character(forwardRepeat[[1]])
    searchRange <- c(forwardRepeat$start + elementLength[1], forwardRepeat$start + elementLength[2])

    if (searchRange[2] > length(Genome[Genome@ranges@NAMES == chr][[1]])) {
      searchRange[2] <- length(Genome[Genome@ranges@NAMES == chr][[1]])
    }

    # consider all reverse matches in range with matching TSD sequences
    reverseRepeats <- reverseMatches[reverseMatches$seqnames == as.character(forwardRepeat$seqnames) &
      reverseMatches$end > searchRange[1] &
      reverseMatches$end < searchRange[2] &
      reverseMatches$strand == "-" &
      reverseMatches$TSD == as.character(forwardRepeat$TSD), ]

    # append matches to packMatches
    if (length(reverseRepeats[, 1]) > 0) {
      for (reverseMatch in 1:length(reverseRepeats[, 1])) {
        packMatches <- rbind(
          packMatches,
          data.frame(
            seqnames = forwardRepeat$seqnames,
            start = forwardRepeat$start,
            end = reverseRepeats[reverseMatch, ]$end,
            width = reverseRepeats[reverseMatch, ]$end - forwardRepeat$start + 1,
            strand = "*"
          )
        )
      }
    }
  }
  return(packMatches)
}

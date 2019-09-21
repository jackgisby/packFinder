#' @title Get Flanking Terminal Site Duplications
#' @description Gets the flanking TSD sequences of the TIRs in \code{tirMatches}.
#' @param tirMatches A dataframe containing genomic ranges and names referring to TIR sequences.
#' @param Genome A DNAStringSet object containing sequences referred to in \code{tirMatches}.
#' @param tsdLength The length of the TSD region to be retrieved (integer).
#' @param strand The strand of the TIR; "+" for forward, "-" for reverse.
#' @author Jack Gisby
#' @details
#' Called by \code{\link{packSearch}}. It is recommended to use the general pipeline function \code{\link{packSearch}} for identification of potential pack elements, which returns TSD sequences as a feature of results, however each stage may be called individually.
#' @return The dataframe \code{tirMatches} with an additional \code{TSD} feature containing flanking TSD sequences as characters.
#' @export

getTsds <- function(tirMatches,
                    Genome,
                    tsdLength,
                    strand) {
  if (strand == "+") {
    tirMatches <- tirMatches[tirMatches$start > tsdLength, ]
    tirMatches$TSD <- mapply(function(seqnames, start, tsdLength, Genome) {
      return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(start - tsdLength):(start - 1)]))
    },
    tirMatches$seqnames,
    tirMatches$start,
    MoreArgs = list(tsdLength = tsdLength, Genome = Genome)
    )

    return(tirMatches)
  } else if (strand == "-") {
    removeMatch <- mapply(function(end, seqnames, tsdLength, Genome) {
      return((end + tsdLength) > Genome[Genome@ranges@NAMES == seqnames][[1]]@length)
    },
    tirMatches$end,
    tirMatches$seqnames,
    MoreArgs = list(tsdLength, Genome)
    )
    tirMatches <- tirMatches[removeMatch == FALSE, ]
    tirMatches$TSD <- mapply(function(seqnames, end, tsdLength, Genome) {
      return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(end + 1):(end + tsdLength)]))
    },
    tirMatches$seqnames,
    tirMatches$end,
    MoreArgs = list(tsdLength = tsdLength, Genome = Genome)
    )
    return(tirMatches)
  }
}

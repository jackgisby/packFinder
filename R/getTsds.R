#' @title Get Flanking Terminal Site Duplications
#' @description Gets the flanking TSD sequences of the TIRs in \code{tirMatches}.
#' @param tirMatches A dataframe containing genomic ranges and names referring to TIR sequences.
#' @param Genome A DNAStringSet object containing sequences referred to in \code{tirMatches}.
#' @param tsdLength The length of the TSD region to be retrieved (integer).
#' @param direction The direction of the TIR; "+" for forward, "-" for reverse.
#' @author Jack Gisby
#' @details
#' Called by \code{\link{packSearch}}. Function intended for internal use.
#' @return The dataframe \code{tirMatches} with an additional \code{TSD} feature containing flanking TSD sequences as characters.


getTSDs <- function(tirMatches,
                    Genome,
                    tsdLength,
                    direction) {

  if(direction == "+") {
    return(tirMatches %>%
             dplyr::filter(start > tsdLength) %>%
             dplyr::mutate(TSD = mapply(function(seqnames, start, tsdLength, Genome) {
               return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(start - tsdLength):(start - 1)]))},
               seqnames,
               start,
               MoreArgs = list(tsdLength = tsdLength, Genome = Genome))
             ))

  } else if (direction == "-") {
    return(tirMatches %>%
             dplyr::mutate(removeMatch = mapply(function(end, seqnames, tsdLength, Genome)  {
               return((end + tsdLength) > Genome[Genome@ranges@NAMES == seqnames][[1]]@length)},
               end,
               seqnames,
               MoreArgs = list(tsdLength, Genome))) %>%
             dplyr::filter(removeMatch == FALSE) %>%
             dplyr::select(-c(removeMatch)) %>%
             dplyr::mutate(TSD = mapply(function(seqnames, end, tsdLength, Genome) {
               return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(end + 1):(end + tsdLength)]))},
               seqnames,
               end,
               MoreArgs = list(tsdLength = tsdLength, Genome = Genome))))
  }
}

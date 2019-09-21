#' @title Extract Sequences from Ranges
#' @description Gets the sequences referred to in \code{packMatches} and returns the sequences added as an additional column to the dataframe..
#' @param packMatches A dataframe containing genomic ranges and names referring to sequences to be extracted. Obtained from \code{\link{packSearch}}.
#' @param Genome A DNAStringSet object containing sequences referred to in \code{packMatches}.
#' @param output The type of object to be returned. At default, returns a \code{\link[Biostrings]{DNAStringSet}} object; alternatively, can return a \code{character} vector.
#' @author Jack Gisby
#' @return The transposon sequences extracted from \code{packMatches}.
#' @export

getPackSeqs <- function(packMatches,
                        Genome,
                        output = "DNAStringSet") {

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

identifyTirMatches <- function(tirMatches, subSeq, Genome, mismatch, strand = "*") {
  # searches genome for potential TIRs based on sequence similarity
  #
  # ---input---
  # tirMatches: dataframe of previously identified matches
  # subSeq: DNAString object containing the sub-sequence to be searched for
  # Genome: DNAStringSet object containing the genome to be searched
  # mismatch: numeric, acceptable edit distance between subSeq and a given sequence in Genome
  # strand: string, direction of match - forward (+) or reverse (-)
  #
  # ---returns---
  # tirMatches: dataframe of previously identified matches and matches identified during this search


  for(i in 1:length(Genome)) {
    matches <- matchPattern(subSeq, Genome[[i]], max.mismatch = mismatch, with.indels = TRUE)

    if(length(matches) > 0) {
      tirMatches <- rbind(tirMatches, data.frame(seqnames = names(Genome)[i],
                                                 start = matches@ranges@start,
                                                 end = matches@ranges@start + matches@ranges@width - 1,
                                                 width = matches@ranges@width,
                                                 strand = strand))
    }
  }

  return(tirMatches)
}

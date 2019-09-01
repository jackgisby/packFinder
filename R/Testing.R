forwardMatches <- as.data.frame(GRanges()) %?%
  identifyTIRMatches(forwardMatches, subSeq, Genome, mismatch = 2, strand = "+") %>%
  getTSDs(forwardMatches, Genome, direction = "+") %>%
  mutate(forwardMatches, TSD = mapply(function(seqnames, start, TSD.length, Genome) {
    return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(end + 1):(end + TSD.length)]))},
    seqnames,
    start,
    MoreArgs = list(TSD.length = TSD.length, Genome = Genome)))


getTSDs <- function(TIR_Matches, Genome, TSD.length, direction) {
  
  if(direction == "+") {
    return(TIR_Matches %>% mutate(
      TSD = mapply(function(seqnames, start, TSD.length, Genome) {
        return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(start - TSD.length):(start - 1)]))},
        seqnames,
        start,
        MoreArgs = list(TSD.length = TSD.length, Genome = Genome))
    ))
  } else if(direction == "-") {
    return(TIR_Matches %>% mutate(forwardMatches, TSD = mapply(function(seqnames, start, TSD.length, Genome) {
      return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(end + 1):(end + TSD.length)]))},
      seqnames,
      start,
      MoreArgs = list(TSD.length = TSD.length, Genome = Genome))))
  }
}
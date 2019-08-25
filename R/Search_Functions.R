identifyTIRMatches <- function(TIR_Matches, subSeq, Genome, mismatch, strand = "*") {
  
  for(i in 1:length(Genome)) {
    matchData <- GRanges(seqnames = names(Genome)[i], 
                         ranges = as(matchPattern(subSeq, Genome[[i]], max.mismatch = mismatch, with.indels = TRUE), "IRanges"),
                         strand = strand)
    
    if(!is.null(matchData)) {
      TIR_Matches <- c(TIR_Matches, matchData)
    }
  }
  
  return(TIR_Matches)
}

identifyPotentialPackElements <- function(subSeq, Genome, element.length, TIR_Matches) {
  potTransposons <- TIR_Matches[strand(TIR_Matches) == "+"]
  reverseMatchess <- TIR_Matches[strand(TIR_Matches) == "-"]
  
  for(forwardMatch in 1:length(potTransposons)) {
    print("yay")
  }
  
  return(potTransposons)
}
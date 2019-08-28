identifyTIRMatches <- function(TIR_Matches, subSeq, Genome, mismatch, strand = "*") {
  
  for(i in 1:length(Genome)) { #refactor without Granges?
    matchData <- GRanges(seqnames = names(Genome)[i], 
                         ranges = as(matchPattern(subSeq, Genome[[i]], max.mismatch = mismatch, with.indels = TRUE), "IRanges"),
                         strand = strand)
    
    if(!is.null(matchData)) {
      TIR_Matches <- rbind(TIR_Matches, as.data.frame(matchData))
    }
  }
  
  return(TIR_Matches)
}

identifyPotentialPackElements <- function(forwardMatches, reverseMatches, subSeq, Genome, mismatch, element.length, TSD.length) {
  
  potTransposons <- as.data.frame(GRanges())
  
  for(forwardMatch in 1:length(forwardMatches[,1])) { #for each forward match
    
    forwardRepeat <- forwardMatches[forwardMatch,]
    chr <- as.character(forwardRepeat[[1]])
    searchRange <- c(forwardRepeat$start + element.length[1], forwardRepeat$start + element.length[2])
    
    if(searchRange[2] > length(Genome[Genome@ranges@NAMES == chr][[1]])) {
      searchRange[2] <- length(Genome[Genome@ranges@NAMES == chr][[1]])
    }
    
    reverseRepeats <- filter(reverseMatches,
                             seqnames == forwardRepeat$seqnames & 
                             end > searchRange[1] & 
                             end < searchRange[2] & 
                             strand == "-"  &
                             TSD == forwardRepeat$TSD)

    if(length(reverseRepeats[,1]) > 0) { #replacing previous loop with df filtering
      for(reverseMatch in 1:length(reverseRepeats[,1])) {
        potTransposons <- rbind(potTransposons, data.frame(seqnames = forwardRepeat$seqnames,
                                                             start = forwardRepeat$start,
                                                             end = reverseRepeats[reverseMatch,]$end,
                                                             width = reverseRepeats[reverseMatch,]$end - forwardRepeat$start,
                                                             strand = "*"))
      }
    }
  }
  return(potTransposons)
}

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
    return(TIR_Matches %>% mutate(
      TSD = mapply(function(seqnames, end, TSD.length, Genome) {
        return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(end + 1):(end + TSD.length)]))},
        seqnames,
        end,
        MoreArgs = list(TSD.length = TSD.length, Genome = Genome))))
  }
}
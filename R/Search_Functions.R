identifyTIRMatches <- function(TIR_Matches, subSeq, Genome, mismatch, strand = "*") {
  
  for(i in 1:length(Genome)) {
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
  
  for(forwardMatch in 1:length(forwardMatches[,1])) {
    
    forwardRepeat <- forwardMatches[forwardMatch,]
    chr <- as.character(forwardRepeat[[1]])
    searchRange <- forwardRepeat$start
    searchRange <- c(searchRange + element.length[1], searchRange + element.length[2])
    
    if(searchRange[2] > length(Genome[Genome@ranges@NAMES == chr][[1]])) {
      searchRange[2] <- length(Genome[Genome@ranges@NAMES == chr][[1]])
    }
    
    reverseRepeats <- reverseMatches[reverseMatches$seqnames == forwardRepeat$seqnames
                                     & reverseMatches$end > searchRange[1]
                                     & reverseMatches$end < searchRange[2]
                                     & reverseMatches$strand == "-",]

    if(length(reverseRepeats[,1]) > 0) {
      fTSD <- (forwardRepeat$start - TSD.length):(forwardRepeat$start - 1) #todo: will cause a bug if flank is out of bounds
      
      for(reverseMatch in 1:length(reverseRepeats[,1])) {
        rTSD <- (reverseRepeats[reverseMatch,]$end + 1):(reverseRepeats[reverseMatch,]$end + TSD.length) #todo: will cause a bug if flank is out of bounds

        if(Genome[Genome@ranges@NAMES == chr][[1]][fTSD] == Genome[Genome@ranges@NAMES == chr][[1]][rTSD]) {
          
          potTransposons <- rbind(potTransposons, data.frame(seqnames = forwardRepeat$seqnames,
                                                             start = forwardRepeat$start,
                                                             end = reverseRepeats[reverseMatch,]$end,
                                                             width = reverseRepeats[reverseMatch,]$end - forwardRepeat$start,
                                                             strand = "*"))
        }
      }
    }
  }
  return(potTransposons)
}
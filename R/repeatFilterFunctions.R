filterElements <- function(potentialPacks, repeatMaps) {
  isTransposon <- vector("list", length(potentialPacks[,1]))
  potentialPacksGRanges <- makeGRangesFromDataFrame(potentialPacks)
  
  for(i in 1:length(isTransposon)) {
    isTransposon[i] <- countOverlaps(potentialPacksGRanges[i], repeatMaps)
    if(countOverlaps(potentialPacksGRanges[i], repeatMaps) > 0) {
      isTransposon[i] <- TRUE
    } else {
      isTransposon[i] <- FALSE
    }
  }
  
  potentialPacks$isTransposon <- isTransposon
  return(potentialPacks)
}

#filter using repeatmap
# repeatMaps <- getRepeatMaps(Genome)
# potentialPacks <- filterElements(potentialPacks, repeatMaps)
# 
# knownCACTA <- saveReport(potentialPacks, subSeq, Genome, integrityFilter = NULL, mismatch = max.mismatch)
# knownCACTA <- saveReport(filter(potentialPacks, isTransposon == FALSE), subSeq, Genome, integrityFilter = NULL, mismatch = max.mismatch)
# print(end-start)
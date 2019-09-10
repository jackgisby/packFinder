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

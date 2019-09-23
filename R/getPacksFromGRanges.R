
getPacksFromGRanges <- function(packGRanges, Genome = NULL, tsdLength = NULL) {
  if (is.null(Genome) | is.null(tsdLength)) {
    return(as.data.frame(packGRanges))
  }
  else if (!is.null(Genome) & !is.null(tsdLength)) {
    return(getTsds(as.data.frame(packGRanges), Genome, tsdLength))
  }
}

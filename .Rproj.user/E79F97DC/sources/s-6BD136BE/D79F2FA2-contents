tirHClust <- function(potentialPacks, Genome, model = "K80", plotDendrogram = TRUE) {
  # gets clusters from the TIRs in potentialPacks and the model specified by user (see ape::dist.dna)
  # requires additional info on the genome of origin
  
  TIRs <- DNAStringSet(c(as.character(potentialPacks$forward_TIR), as.character(potentialPacks$reverse_TIR)))
  TIRs@ranges@NAMES <- c(paste0(potentialPacks$Genome, "f", row.names(potentialPacks)), paste0(potentialPacks$Genome, "r", row.names(potentialPacks)))

  clust <- as.DNAbin(TIRs) %>%
    #dist.dna(model = model)
    kdistance(k=7)
  clust[is.na(clust)] <- 0
  clust[is.nan(clust)] <- 0
  clust[is.infinite(clust)] <- 0.75
  
  dend <- clust %>%
    hclust() %>%
    as.dendrogram() %>%
    highlight_branches_col()
  
  dirCols <- ifelse(grepl("f", labels(dend)), 3, 4)
  labels_colors(dend) <- dirCols
  
  if(plotDendrogram == TRUE) {
    png("Data/Output/Plots/TIR_Relationships.png", width = 1000, height = 500)
    plot(dend, main = "TIR Relationships")
    colored_bars(colors = dirCols, dend=dend, sort_by_labels_order = FALSE, rowLabels = "direction")
    dev.off()
    
    plot(dend, main = "TIR Relationships")
    colored_bars(colors = dirCols, dend=dend, sort_by_labels_order = FALSE, rowLabels = "direction")
  }
  return(dend)
}

getTirConsensus <- function(potentialPacks, dend, h) {
  # takes dendrogram and separates into clusters to produce consensus sequences of TIRs
  clust <- cutree(dend, h = h)
  consensusSeqs <- vector("list", length = length(unique(clust)))
  
  for(i in 1:length(unique(clust))) {
    seqNames <- names(clust[clust == unique(clust)[i]])
    dir <- grepl("f", seqNames) 
    
    consensusSeqs[[i]] <-
      c(DNAStringSet(potentialPacks$forward_TIR[as.integer(subseq(seqNames[dir==TRUE], start = 3))]),
      DNAStringSet(potentialPacks$reverse_TIR[as.integer(subseq(seqNames[dir==FALSE], start = 3))])) %>%
      consensusString()
    
    names(consensusSeqs)[i] <- paste(seqNames, collapse = ', ')
  }
  
  return(consensusSeqs)
}

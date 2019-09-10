getTirClusters <- function(potentialPacks, model = "TN93") {
  # gets clusters from the TIRs in potentialPacks and the model specified by user (see ape::dist.dna)
  
  TIRs <- DNAStringSet(c(potentialPacks$forward_TIR, potentialPacks$reverse_TIR))
  TIRs@ranges@NAMES <- c(paste0("f", row.names(potentialPacks)), paste0("r", row.names(potentialPacks)))
  
  dend <- as.DNAbin(TIRs) %>%
    dist.dna(model = model) %>% 
    njs() %>%
    as.dendrogram() %>%
    highlight_branches_col() %>%
    hang.dendrogram(hang = -1) %>%
    sort(type="nodes")
  
  #idCols <- mapply(function(x) {return(as.integer(substr(x, 2, 3)))}, labels(dend), SIMPLIFY = TRUE)
  dirCols <- ifelse(grepl("f", labels(dend)), 3, 4)
  labels_colors(dend) <- dirCols
  
  png("Data/Output/Plots/TIR_Relationships.png", width = 1000, height = 500)
  plot(dend, main = "TIR Relationships")
  colored_bars(colors = dirCols, dend=dend, sort_by_labels_order = FALSE, rowLabels = "direction")
  dev.off()
  
  plot(dend, main = "TIR Relationships")
  colored_bars(colors = dirCols, dend=dend, sort_by_labels_order = FALSE, rowLabels = "direction")
  
  return(dend)
}

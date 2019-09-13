getTirHClust <- function(potentialPacks, model = "K80") {
  # gets clusters from the TIRs in potentialPacks and the model specified by user (see ape::dist.dna)
  
  TIRs <- DNAStringSet(c(as.character(potentialPacks$forward_TIR), as.character(potentialPacks$reverse_TIR)))
  TIRs@ranges@NAMES <- c(paste0("f", row.names(potentialPacks)), paste0("r", row.names(potentialPacks)))
  
  dend <- as.DNAbin(TIRs) %>%
    dist.dna(model = model) %>% 
    hclust() %>%
    as.dendrogram() %>%
    highlight_branches_col()

  #idCols <- mapply(function(x) {return(as.integer(substr(x, 2, 3)))}, labels(dend), SIMPLIFY = TRUE)
  dirCols <- ifelse(grepl("f", labels(dend)), 3, 4)
  labels_colors(dend) <- dirCols
  organismCols <- mapply(function(name, uniqueNames) {
    return(which(uniqueNames == name) + 2)},
    organism,
    MoreArgs = list(unique(organism)))
  

  png("Data/Output/Plots/TIR_Relationships.png", width = 1000, height = 500)
  plot(dend, main = "TIR Relationships")
  colored_bars(colors = dirCols, dend=dend, sort_by_labels_order = FALSE, rowLabels = "direction")
  dev.off()

  plot(dend, main = "TIR Relationships")
  colored_bars(colors = cbind(dirCols, organismCols), dend=dend, sort_by_labels_order = FALSE, rowLabels = c("direction", "organism"))

  return(cbind(dirCols, organismCols))
}

getOrganismHClust <- function(potentialPacks, genomeList, model = "K80") {
  # gets clusters from the TIRs in potentialPacks and the model specified by user (see ape::dist.dna)
  # requires additional info on the genome of origin
  
  potentialPacks$Genome <- mapply(function(name, uniqueNames) {
    return(which(uniqueNames == name) + 4)},
    genomeList,
    MoreArgs = list(unique(genomeList)))
  TIRs <- DNAStringSet(c(as.character(potentialPacks$forward_TIR), as.character(potentialPacks$reverse_TIR)))
  TIRs@ranges@NAMES <- c(paste0(potentialPacks$Genome, "f", row.names(potentialPacks)), paste0(potentialPacks$Genome, "r", row.names(potentialPacks)))

  clust <- as.DNAbin(TIRs) %>%
    #dist.dna(model = model)
    kdistance(k=7)
  clust[is.na(clust)] <- 0
  clust[is.nan(clust)] <- 0
  clust[is.infinite(clust)] <- 1
  
  dend <- clust %>%
    hclust() %>%
    as.dendrogram() %>%
    highlight_branches_col()
  
  #idCols <- mapply(function(x) {return(as.integer(substr(x, 2, 3)))}, labels(dend), SIMPLIFY = TRUE)
  dirCols <- ifelse(grepl("f", labels(dend)), 3, 4)
  labels_colors(dend) <- dirCols
  organismCols <- lapply(labels(dend), function(x) subseq(x, start=1, end=1))
  
  png("Data/Output/Plots/TIR_Relationships.png", width = 1000, height = 500)
  plot(dend, main = "TIR Relationships")
  colored_bars(colors = dirCols, dend=dend, sort_by_labels_order = FALSE, rowLabels = "direction")
  dev.off()
  
  plot(dend, main = "TIR Relationships")
  colored_bars(colors = organismCols, dend=dend, sort_by_labels_order = FALSE, rowLabels = "organism")
  legend("topright", legend = unique(genomeList), fill = 5:(length(unique(genomeList))+5))
  
  return(clust)
}

getOrganismKClust <- function(potentialPacks, genomeList, model = "K80") {
  # gets clusters from the TIRs in potentialPacks and the model specified by user (see ape::dist.dna)
  # requires additional info on the genome of origin
  
  potentialPacks$Genome <- mapply(function(name, uniqueNames) {
    return(which(uniqueNames == name) + 4)},
    genomeList,
    MoreArgs = list(unique(genomeList)))
  TIRs <- DNAStringSet(c(as.character(potentialPacks$forward_TIR), as.character(potentialPacks$reverse_TIR)))
  TIRs@ranges@NAMES <- c(paste0("f", row.names(potentialPacks)), paste0("r", row.names(potentialPacks)))
  
  clust <- as.DNAbin(TIRs) %>%
    #dist.dna(model = model)
    kdistance()
  
  cmdscale(clust, eig=TRUE) %>%
    fortify() %>%
    setNames(c("MDS1", "MDS2")) %>%
    ggplot(aes(MDS1, MDS2, colour = as.character(rep(potentialPacks$Genome, 2)),
               label = TIRs@ranges@NAMES)) +
    geom_point() + 
    # geom_text(aes(label = TIRs@ranges@NAMES), position = position_nudge(x = 0.02), colour = "black") +
    # scale_color_manual(labels = "black") +
    scale_colour_discrete(name = "Organism", labels = unlist(unique(genomeList))) %>%
    return()
}



  

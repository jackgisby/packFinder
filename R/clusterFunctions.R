packClust <- function(packMatches, 
                        Genome, 
                        identity = 0.6, 
                        threads = 1, 
                        strand = "both", 
                        saveFolder = "packFinderData/Output/vSearch/", 
                        vSearchPath = "D:/vsearch-2.14.1-win-x86_64/vsearch.exe") {
  
  
  packMatchesFile <- "packFinderData/Data/packMatches.fa"
  packMatches <- packMatches %>% 
    mutate(ID = 1:length(packMatches[,1])) %>%
    arrange(desc(width))
  
  packMatchesSet <- DNAStringSet(packMatches$seq)
  packMatchesSet@ranges@NAMES <- as.character(packMatches$ID)
  writeXStringSet(packMatchesSet, packMatchesFile)
  
  system2(
    command = vSearchPath,
    args = paste0(
      "--cluster_smallmem ",
      packMatchesFile,
      " \ ",
      "--qmask none \ ",
      "--uc ",
      file.path(saveFolder, paste0("packMatches", ".uc")),
      " \ ",
      "--id ",
      identity,
      " \ ",
      "--threads ",
      threads,
      " \ ",
      "--clusterout_sort \ ",
      "--clusterout_id \ ",
      "--strand ",
      strand,
      " \ ",
      "--log ",
      file.path(saveFolder, paste0("packMatches", ".log")),
      " \ ",
      "--blast6out ",
      file.path(saveFolder, paste0("packMatches", ".blast6out")),
      " \ ",
      "--sizeout"
    )
  )
  
  vSearchClusts <- read.uc(file.path(saveFolder, paste0("packMatches", ".uc"))) %>%
    filter(type != "C") %>%
    arrange(query)
  packMatches %>%
    mutate(strand = mapply(function(strand) {
      if(strand == "*") {
        return("+")
      } else {
        return(strand)
      }
    },
    strand = as.character(vSearchClusts$strand))) %>%
    mutate(cluster = vSearchClusts$cluster) %>%
      return()
}

read.uc <- function(savePath) {
  packClusts <- read.table(savePath, sep = "\t") 
  colnames(packClusts) <- c("type",
                            "cluster",
                            "width",
                            "identity",
                            "strand",
                            "6",
                            "7",
                            "cigarAlignment",
                            "query",
                            "target"
  )
  
  packClusts %>%
    select(-c("6", "7")) %>%
    return()
}

tirClust <- function(packMatches, plot = TRUE, plotSavePath = NULL, k = 5, tirLength = 25) {
  fConsensusSeqs <- vector("list", length = length(unique(packMatches$clustID)))
  rConsensusSeqs <- vector("list", length = length(unique(packMatches$clustID)))
  
  for(c in 1:length(unique(packMatches$cluster))) {
    clustID <- unique(packMatches$cluster)[c]
    clust <- filter(packMatches, cluster == clustID)
    forwardTirs <- vector("list", length = length(clust[,1]))
    reverseTirs <- vector("list", length = length(clust[,1]))

    for(i in 1:length(clust[,1])) {
      if(clust$strand[i] == "+") {
        forwardTirs[[i]] <- Genome[Genome@ranges@NAMES == clust$seqnames[i]][[1]][clust$start[i]:(clust$start[i]+tirLength)]
        reverseTirs[[i]] <- reverseComplement(Genome[Genome@ranges@NAMES == clust$seqnames[i]][[1]][(clust$end[i]-tirLength):clust$end[i]])
      } else if(clust$strand[i] == "-") {
        forwardTirs[[i]] <- reverseComplement(Genome[Genome@ranges@NAMES == clust$seqnames[i]][[1]][(clust$end[i]-tirLength):clust$end[i]])
        reverseTirs[[i]] <- Genome[Genome@ranges@NAMES == clust$seqnames[i]][[1]][clust$start[i]:(clust$start[i]+tirLength)]
      }
    }

    fConsensusSeqs[[c]] <- consensusString((DNAStringSet(forwardTirs)))
    rConsensusSeqs[[c]] <- consensusString((DNAStringSet(reverseTirs)))
  }
  
  fConsensusSeqs <- DNAStringSet(unlist(fConsensusSeqs))
  rConsensusSeqs <- DNAStringSet(unlist(rConsensusSeqs))
  fConsensusSeqs@ranges@NAMES <- paste0("f", as.character(unique(packMatches$cluster)))
  rConsensusSeqs@ranges@NAMES <- paste0("r", as.character(unique(packMatches$cluster)))
  consensusSeqs <- c(fConsensusSeqs, rConsensusSeqs)
  
  if(plot == TRUE) {
    dend <- as.DNAbin(consensusSeqs) %>%
      kdistance(k=k) %>%
      hclust() %>%
      as.dendrogram() %>%
      highlight_branches_col()
    
    plot(dend, main = "Clustered Transposon TIR Relationships")
  }
  
  if(!is.null(plotSavePath)) {
    png(plotSavePath, width = 1000, height = 600)
    plot(dend, main = "Clustered Transposon TIR Relationships")
    dev.off()
  }
  
  return(consensusSeqs)
}

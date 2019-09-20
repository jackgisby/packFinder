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

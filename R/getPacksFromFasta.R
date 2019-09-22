getPacksFromFasta <- function(file) {
  fileCon <- file(file, "r")
  packMatches <- data.frame(
    seqnames = character(),
    start = integer(),
    end = integer(),
    width = integer(),
    strand = character(),
    TSD = character()
  )
  while (TRUE) {
    seqName <- readLines(fileCon, n = 1)
    if ((length(seqName) == 0) | (length(seq) == 0)) {
      break
    } else if (substr(seqName, 1, 1) != ">") {
      break
    }
    seq <- readLines(fileCon, n = 1)

    seqName <- gsub(">","",seqName)
    seqName <- gsub("start =","",seqName)
    seqName <- gsub("end =","",seqName)
    seqName <- gsub("width =","",seqName)
    seqName <- gsub("strand =","",seqName)
    seqName <- gsub("TSD =","",seqName)
    seqName <- gsub(" ","",seqName)
    seqName <- strsplit(seqName, "|", fixed = TRUE)[[1]]
    names(seqName) <- c("seqnames", "start", "end", "width", "strand", "TSD")
    packMatches <- rbind(packMatches,
                         seqName,
                         stringsAsFactors = FALSE)
    seqName <- NULL
  }
  close(fileCon)
  colnames(packMatches) <- c("seqnames", "start", "end", "width", "strand", "TSD")
  return(packMatches)
}

getPacksFromFasta("devData/test.fasta")

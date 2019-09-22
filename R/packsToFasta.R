packsToFasta <- function(packMatches, file, Genome) {
  file.create(file)
  for(match in 1:length(packMatches[,1])) {
    write(c(paste0(">", packMatches$seqnames[match],
                       " | start = ", packMatches$start[match],
                       " | end = ", packMatches$end[match],
                       " | width = ", packMatches$width[match],
                       " | strand = ", packMatches$strand[match],
                       " | TSD = ", packMatches$TSD[match]),
                 as.character(Genome[Genome@ranges@NAMES == packMatches$seqnames[match]][[1]][packMatches$start[match]:packMatches$end[match]])),
          file = file,
          append = TRUE)
  }
  return(print(paste0("FASTA written to ", file)))
}

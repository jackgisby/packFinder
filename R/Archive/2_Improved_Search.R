SubSeq = strsplit("CACTACAAGAAATATGACATTGGTA", "")[[1]]
setwd("C:/Users/jackg/OneDrive/Documents/Work/Final Year 2019-2020/Transposon_Search")

SubSeq_Length = length(SubSeq)
SearchSeq = strsplit(readChar("Example_Seq.txt", file.info("Example_Seq.txt")$size), "")[[1]]
SearchSeq_Length = length(SearchSeq)

j = 1
i = 1

start = Sys.time()

for (i in 1:(SearchSeq_Length - SubSeq_Length)) {
  
  if (SearchSeq[i] == SubSeq[j]) { 
    j = j + 1
    if (j == SubSeq_Length + 1) {
      print(paste("i = ", i-24))
      print(paste("p = ", ((1/(4^SubSeq_Length))*SearchSeq_Length)))
      j = 1
    }
    
  } else {
    if (SearchSeq[i] == SubSeq[1]) {
      j = 2
    } else {
      j = 1
    }
  }
}

end = Sys.time()

end - start
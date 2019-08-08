SubSeq = DNAString("CACTACAAGAAATATGACATTGGTA")
setwd("C:/Users/jackg/OneDrive/Documents/Work/Final Year 2019-2020/Transposon_Search")
file = "Example_Seq.txt"
SearchSeq = readChar(file, file.info(file)$size)
SearchSeq = b
start = Sys.time()
for (i in 1:(nchar(SearchSeq) - nchar(SubSeq))) {
  if (substr(SearchSeq, i, i+nchar(SubSeq)-1) == SubSeq) {
    print("found", str(i))
    print(substr(SearchSeq, i, i+nchar(SubSeq)-1))
  }
}


end = Sys.time()

end - start
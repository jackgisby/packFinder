library(Biostrings)
library(BSgenome.Athaliana.TAIR.TAIR9)

source("R/Pipeline.R")

SubSeq <- DNAString("CACTACAAAAATATCATTTTA")
SubSeq <- DNAString("CTAGTCATG")
ArAth <- BSgenome.Athaliana.TAIR.TAIR9

start = Sys.time()
Search_Visualise(SubSeq, ArAth, max.mismatch = 0, with.indels = FALSE, fixed = TRUE)
end = Sys.time()
print(end-start)
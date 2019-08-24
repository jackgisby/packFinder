library(Biostrings)
library(BSgenome.Athaliana.TAIR.TAIR9)

source("R/Pipeline.R")

SubSeq <- "CACTACAAAAATATCATTTTA"
ArAth <- BSgenome.Athaliana.TAIR.TAIR9

start = Sys.time()
#Search_Visualise(SubSeq, ArAth, max.mismatch = 2, with.indels = FALSE, fixed = TRUE)
end = Sys.time()
print(end-start)
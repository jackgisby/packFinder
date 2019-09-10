library(Biostrings)
library(biomartr)
library(GenomicRanges)
library(dplyr)
library(rBLAST)
library(hoardeR)

for(file in 1:length(list.files("R/"))) {
  source(paste0("R/", list.files("R/")[file]))
}
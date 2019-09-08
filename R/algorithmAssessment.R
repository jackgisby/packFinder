source("R/packSearch.R")
source("R/devFunctions.R")

#pack find
Genome <- initialise()
subSeq <- DNAString("CACTACAA") #CACTACAA-AAATAT / DNAString(consensusString(knownTIRs))
max.mismatch = 0

start <- Sys.time()
potentialPacks <- packSearch(subSeq, Genome, mismatch = max.mismatch, element.length = c(300, 3500), TSD.length = 3)
end <- Sys.time()

knownCACTA <- saveReport(potentialPacks, subSeq, Genome, integrityFilter = NULL, mismatch = max.mismatch)
print(end-start)

db <- blast(db="C:/Users/jackg/Documents/R/nt_db/nt/nt", type = "blastn")

#blast
start <- Sys.time()
#BlastMatches <- packBlast(potentialPacks, "nt")
#blastMatches <- blastSeq(DNAStringSet(Genome$`NC_003070.9 Arabidopsis thaliana chromosome 1 sequence`[13131404:13132121]),database = "nt")
blastMatch <- predict(db, 
                      DNAStringSet(Genome$`NC_003070.9 Arabidopsis thaliana chromosome 1 sequence`[13131404:13132121]), 
                      BLAST_args = "-num_threads 4")
end <- Sys.time()
print(end-start)


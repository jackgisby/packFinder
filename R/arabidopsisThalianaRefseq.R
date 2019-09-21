#' @title Arabidopsis thaliana Refseq Genome
#' @format A DNAStringSet object of length 7.
#' @docType data
#' @usage data(arabidopsisThalianaRefseq)
#' @example \dontrun{
#' data(arabidopsisThalianaRefseq)
#'
#' subSeq <- Biostrings::DNAString("CACTACAA")
#' packMatches <- packSearch(subSeq, arabidopsisThalianaRefseq, elementLength = c(300,3500), tsdLength = 3)
#' }
#' @source The Arabidopsis thaliana genome was downloaded from the NCBI refseq database on 20/SEP/2019 using \code{\link[biomartr]{getGenome}}. The genome can be accessed from the NCBI ftp server: \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz}.

"arabidopsisThalianaRefseq"

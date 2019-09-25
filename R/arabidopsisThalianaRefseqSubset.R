#' @title
#' Arabidopsis thaliana Refseq Genome Chromosome 3
#'
#' @description
#' The chromosome 3 reference sequence for Arabidopsis thaliana as a
#' \code{\link[Biostrings]{DNAStringSet}}. Can be used as a test data set,
#' as in the associated introduction vignette.
#'
#' @format
#' A DNAStringSet object containing a DNAString for Arabidopsis thaliana's
#' chromosome 3 sequence.
#'
#' @usage
#' data(arabidopsisThalianaRefseq)
#'
#' @examples \dontrun{
#' data(arabidopsisThalianaRefseq)
#'
#' packMatches <- packSearch(
#'   Biostrings::DNAString("CACTACAA"),
#'   arabidopsisThalianaRefseq,
#'   elementLength = c(300, 3500),
#'   tsdLength = 3)}
#'
#' @source
#' The Arabidopsis thaliana genome was downloaded from the NCBI refseq
#' database on 20/SEP/2019, using \code{\link[biomartr]{getGenome}},
#' and chromosome 3 was extracted. The genome may also be accessed from the
#' NCBI ftp server: \url{ftp://ftp.ncbi.nlm.nih.gov/genomes}.

"arabidopsisThalianaRefseq"

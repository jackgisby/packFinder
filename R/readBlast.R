#' @title
#' Convert NCBI BLAST+ Files to Dataframe
#'
#' @description
#' Reads .blast6out files (NCBI Blast Format) generated by 
#' the VSEARCH clustering and alignment algorithms.
#'
#' @param file
#' The file path of the blast file.
#' 
#' @param minE
#' Blast results with e values greater than
#' the specified cutoff will be ignored.
#' 
#' @param length
#' Blast results alignment lengths lower below
#' this value will be ignored
#' 
#' @param identity
#' Blast results with target sequence identities below
#' this value will be ignored.
#' 
#' @param removeExactMatches
#' If true, matches with 100% sequence identity will 
#' be ignored to prevent self-hits. 
#' 
#' @param scope
#' If specified, blast results below the specified value
#' will be ignored. Note that the dataframe of transposon
#' matches must also be supplied to calculate scope. Scope is 
#' the proportion of the transposon's internal sequence 
#' occupied by the BLAST hit. 
#' 
#' @param packMatches
#' taframe containing genomic ranges and names referring 
#' to sequences to be extracted. Can be obtained from 
#' \code{\link{packSearch}} or generated from a 
#' \code{\link[GenomicRanges:GRanges-class]{GRanges}} object, 
#' after conversion to a dataframe. Must contain the 
#' following features:
#' \itemize{
#'     \item start - the predicted element's start base 
#'     sequence position.
#'     \item end - the predicted element's end 
#'     base sequence position.
#'     \item seqnames - character string 
#'     referring to the sequence name in \code{Genome} to 
#'     which \code{start} and \code{end} refer to.
#' }
#'
#' @return
#' A dataframe containing the converted .blast6out file. 
#' The file contains the following features:
#' \itemize{
#'     \item Query sequence ID
#'     \item Target sequence ID
#'     \item Percenty sequence identity
#'     \item Alignment length
#'     \item Number of mismatches
#'     \item Number of gaps
#'     \item Base position of alignment start 
#'     in query sequence
#'     \item Base position of alignment end in query sequence
#'     \item Base position of alignment start in target sequence
#'     \item Base position of alignment end in target sequence
#'     \item E-value
#'     \item Bit score
#' }
#'
#' @details
#' blast6out file is tab-separated text file compatible with 
#' NCBI BLAST m8 and NCBI BLAST+ outfmt 6 formats. One 
#' cluster/alignment can be found for each line.
#'
#' @seealso
#' code{\link{blastAnalysis}}, code{\link{blastAnnotate}},
#' code{\link{packAlign}}, 
#' code{\link{readUc}}, code{\link{packClust}}
#'
#' @references
#' For further information, see the NCBI BLAST+ application
#' documentation and help pages 
#' (https://www.ncbi.nlm.nih.gov/pubmed/20003500?dopt=Citation).
#' VSEARCH may be downloaded from 
#' \url{https://github.com/torognes/vsearch}; see 
#' \url{https://www.ncbi.nlm.nih.gov/pubmed/27781170} 
#' for further information.
#'
#' @examples 
#' readBlast(system.file(
#'     "extdata", 
#'     "packMatches.blast6out", 
#'     package = "packFinder"
#' ))
#' 
#' @author
#' Jack Gisby
#'
#' @export

readBlast <- function(file, minE = 1, length = 0, identity = 0,
                      removeExactMatches = FALSE, 
                      scope = NULL, packMatches = NULL) {
    
    checkPermissions(file)
    
    blastData <- utils::read.table(file, sep = "\t", 
                           header = FALSE, stringsAsFactors = FALSE)
    
    colnames(blastData) <- c("query_id", "subject_id","identity", 
                             "alignment_length", "mismatches", "gap_opens", 
                             "q.start", "q.end", "s.start", "s.end", "evalue", 
                             "bit_score")
    
    blastData <- blastData[blastData$evalue < minE,]
    blastData <- blastData[blastData$alignment_length > length,]
    blastData <- blastData[blastData$identity > identity,]
    
    if (removeExactMatches) {
        blastData <- blastData[blastData$identity < 100,]
    }
    
    if (!is.null(scope)) {
        width <- vector("integer", nrow(blastData))
        scope_filter <- vector("numeric", nrow(blastData))
        
        for (rowname in seq_len(nrow(blastData))) {
            width[rowname] <- 
                packMatches$width[
                    packMatches$id == blastData[rowname,]$query_id]
            scope_filter[rowname] <- 
                blastData[rowname,]$alignment_length / width[rowname]
        }
        
        blastData$width <- width
        blastData$scope_filter <- scope_filter
        
        blastData <- blastData[blastData$scope_filter > scope,]
    }
    
    return(blastData)
}

#' @title
#' Functional Annotation of PackTYPE Elements
#'
#' @description
#' Uses hits, previously generated using blast, 
#' to annotate transposon hits. Transposons with 
#' non-redundant transposase hits are classed as 
#' autonomous ("auto"), while others are classed as
#' "other" or "pack" based on whether the element 
#' has non-redundant hits to other proteins. 
#'  
#' @param packMatches
#' A dataframe of potential Pack-TYPE transposable elements, 
#' in the format given by \code{\link{packSearch}}. This 
#' dataframe is in the format produced by coercing a 
#' \code{link[GenomicRanges:GRanges-class]{GRanges}} 
#' object to a dataframe: \code{data.frame(GRanges)}. 
#' Will be saved as a FASTA file for VSEARCH.
#' 
#' @param protHits
#' BLAST results for non-transposon related genes or 
#' proteins (as a \code{data.frame}). Generated using 
#' \code{\link{blastAnalysis}}.
#' 
#' @param autoHits
#' BLAST results for transposon related genes or 
#' proteins (as a \code{data.frame}). Generated using 
#' \code{\link{blastAnalysis}}.
#' 
#' @return 
#' Returns the original \code{packMatches} dataframe, 
#' with the addition of a "classification" column 
#' containing one of the following values:
#' \itemize{
#'     \item auto - elements that match known 
#'     transposases or transposon-related proteins 
#'     are classified as autonomous elements
#'     \item pack - elements that match other 
#'     proteins or genic sequences may be classified 
#'     as Pack-TYPE elements
#'     \item other - elements that generate no 
#'     significant hits
#' }
#'     
#' @seealso 
#' \code{\link{blastAnalysis}},
#' \code{\link{readBlast}}, \code{\link{packBlast}}
#' 
#' @examples
#' data("packMatches")
#' 
#' # read in some protein hits
#' p <- data.frame(
#'     query_id = c(2, 3),
#'     subject_id = c("prot", "hyp")
#' )
#' 
#' # read in some autonomous hits
#' a <- data.frame(
#'     query_id = c(3, 4),
#'     subject_id = c("transposase", "mutator")
#' )
#' 
#' blastAnnotate(p, a, packMatches)
#' 
#' @references 
#' For further information, see the NCBI BLAST+ application
#' documentation and help pages 
#' (https://www.ncbi.nlm.nih.gov/pubmed/20003500?dopt=Citation).
#' 
#' @note 
#' Requires that the query ids in the protein and 
#' autonomous hits match the row names in packMatches.
#'  
#' @author
#' Jack Gisby
#'
#' @export

blastAnnotate <- function(protHits, autoHits, packMatches) {
    matchType <- vector("character", length=nrow(packMatches))
    
    for (rowname in seq_len(nrow(packMatches))) {
        protHit <- nrow(protHits[protHits$query_id == as.integer(rowname),])
        autoHit <- nrow(autoHits[autoHits$query_id == as.integer(rowname),])
        
        if (autoHit > 0) {
            matchType[rowname] <- "auto"
        } else if (protHit > 0) {
            matchType[rowname] <- "pack"
        } else {
            matchType[rowname] <- "other"
        }
    }
    
    packMatches$classification <- matchType
    return(packMatches)
}

#' @title
#' Functional Annotation of PackTYPE Elements
#'
#' @description
#'  
#'
#' @param 
#' 
#'     
#' @seealso 
#' \code{\link{packSearch}}
#' 
#' @examples
#' \dontrun{
#' 
#' }
#' 
#' @references 
#' For further information, see the NCBI BLAST+ application
#' documentation and help pages 
#' (https://www.ncbi.nlm.nih.gov/pubmed/20003500?dopt=Citation).
#'  
#' @author
#' Jack Gisby
#'
#' @export

readBlast <- function(hits, cutoff = 1, length = 0, identity = 0,
                      removeExactMatches = FALSE, 
                      Scope = NULL, packMatches = NULL) {
    
    blastData <- read.table(hits, sep = "\t", 
                           header = FALSE, stringsAsFactors = FALSE)
    
    colnames(blastData) <- c("query_id", "subject_id","identity", 
                             "alignment_length", "mismatches", "gap_opens", 
                             "q.start", "q.end", "s.start", "s.end", "evalue", 
                             "bit_score")
    
    blastData <- blastData[blastData$evalue < cutoff,]
    blastData <- blastData[blastData$alignment_length > length,]
    blastData <- blastData[blastData$identity > identity,]
    
    if (removePerfect) {
        blastData <- blastData[blastData$identity < 100,]
    }
    
    if (!is.null(Scope)) {
        width <- vector("integer", nrow(blastData))
        scope <- vector("numeric", nrow(blastData))
        
        for (rowname in 1:nrow(blastData)) {
            width[rowname] <- clusteredMatches$width[clusteredMatches$id == blastData[rowname,]$query_id]
            scope[rowname] <- blastData[rowname,]$alignment_length / width[rowname]
        }
        
        blastData$width <- width
        blastData$scope <- scope
        
        blastData <- blastData[blastData$scope > Scope,]
    }
    
    return(blastData)
}

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

blastAnnotate <- function(protHits, autoHits, packMatches) {
    matchType <- vector("character", length=nrow(packMatches))
    
    for (rowname in 1:nrow(packMatches)) {
        protHit <- nrow(protHits[blastData$query_id == as.integer(rowname),])
        autoHit <- nrow(autoHits[muleHits$query_id == as.integer(rowname),])
        
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

#' @title
#' Pipeline for BLAST/Classification of PackTYPE Elements
#'
#' @description
#' Run BLAST against user-specified databases of 
#' non-transposon and transposon-relates proteins.
#' Can be used to classify transposons based on 
#' their internal sequences. 
#'
#' @param 
#' 
#'     
#' @seealso 
#' \code{\link{packSearch}}
#' 
#' @examples
#' \dontrun{
#' packMatches <- data(packMatches)
#' Genome <- data(arabidopsisThalianaRefseq)
#' packBlast(packMatches, Genome, 
#'     protDb = ", 
#'     autoDb, 
#'     blastPath)
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

packBlast <- function(packMatches, Genome, blastPath, protDb, autoDb,
                      minE = 1e-3, blastTask = "blastn-short", maxHits = 100,
                      threads = 1, saveFolder = NULL, tirCutoff = 100,
                      autoCutoff = 1e-5, autoLength = 150, autoIdentity = 70,
                      autoScope = NULL, protCutoff = 1e-5, protLength = 250, 
                      protIdentity = 70, protScope = 0.3) {
    
    if (is.null(saveFolder)) {
        saveFolder <- getwd()
    }
    
    blastAnalysis(packMatches, Genome, protDb = protDb, autoDb = autoDb, 
                  blastPath = blastPath, minE = minE, blastTask = blastTask,
                  maxHits = maxHits, threads = threads, saveFolder = saveFolder,
                  tirCutoff = tirCutoff)
    
    packMatches$id <- rownames(packMatches)
    
    autoHits <- readBlast(file.path(saveFolder, "autoHits.blast"), 
                          cutoff = autoCutoff, length = autoLength, 
                          identity = autoIdentity, removeExactMatches = FALSE, 
                          Scope = autoScope, packMatches = packMatches)
    
    protHits <- readBlast(file.path(saveFolder, "protHits.blast"), 
                          cutoff = protCutoff, length = protLength, 
                          identity = protIdentity, removeExactMatches = TRUE,
                          Scope = protScope, packMatches = packMatches)
    
    return(blastAnnotate(protHits, autoHits, packMatches))
}

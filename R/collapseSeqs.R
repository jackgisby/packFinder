#' @title
#' Collapse Overlapping Sequences 
#'
#' @description
#' The sequences predicted by \code{\link{packSearch}} often 
#' overlap, which may be due to the presence of closely 
#' interspersed elements or false TIR identification. 
#' In such cases, these elements can be combined using 
#' \code{link[GenomicRanges:GRanges-class]{GRanges}} 
#' in order to collapse overlapping elements, preventing 
#' over-estimation of transposon numbers. Also removes 
#' duplicate elements that have been generated in the 
#' case of multiple searches. 
#'
#' @param packMatches
#' A dataframe containing genomic ranges and names referring 
#' to sequences to be extracted. This dataframe is in the format 
#' produced by coercing a 
#' \code{link[GenomicRanges:GRanges-class]{GRanges}} 
#' object to a dataframe: \code{data.frame(GRanges)}. 
#' 
#' Must contain the following features:
#' \itemize{
#'     \item start - the predicted element's start base 
#'     sequence position.
#'     \item end - the predicted element's end base 
#'     sequence position.
#'     \item seqnames - character string referring to the 
#'     sequence name in \code{Genome} to which \code{start} 
#'     and \code{end} refer to.
#' }
#'
#' @param Genome
#' A DNAStringSet object containing sequences referred to 
#' in \code{packMatches} (the object originally used to 
#' predict the transposons \code{\link{packSearch}}).
#' 
#' @return
#' A set of non-overlapping transposon sequences in the format
#' of the input dataframe. 
#'     
#' @seealso 
#' \code{\link{packSearch}},
#' \code{link[GenomicRanges:GRanges-class]{GRanges}}
#' 
#' @examples
#' data(packMatches)
#' data(arabidopsisThalianaRefseq)
#' 
#' packMatches$start <- 1
#' packMatches$end <- 10
#' 
#' collapseSeqs(packMatches, arabidopsisThalianaRefseq)
#' 
#' @author
#' Jack Gisby
#'
#' @export

collapseSeqs <- function(packMatches, Genome) {
    uniqueMatches <- unique(packMatches[,1:5])
    collapsedMatches <- uniqueMatches[0,]
    
    uniqueRanges <- packsToGRanges(uniqueMatches)
    collapsedRanges <- packsToGRanges(collapsedMatches)
    
    while (length(uniqueRanges) > 0) {
        query <- uniqueRanges[1]
        uniqueRanges[1] <- NULL
        
        overlaps <- GenomicRanges::findOverlaps(query, uniqueRanges, 
                                                ignore.strand = TRUE)
        
        if (length(overlaps) == 0) {
            collapsedRanges <- c(collapsedRanges, query)
        } else {
            queryStart <- start(query)
            queryEnd <- queryStart + GenomicRanges::width(query) - 1
            overlapStart <- start(uniqueRanges[S4Vectors::to(overlaps)[1]])
            overlap <- uniqueRanges[S4Vectors::to(overlaps)[1]]
            overlapEnd <- overlapStart + GenomicRanges::width(overlap) - 1
            
            IRanges::ranges(uniqueRanges[S4Vectors::to(overlaps)[1]]) <- 
                IRanges::IRanges(start = min(queryStart, overlapStart), 
                        end = max(queryEnd, overlapEnd))
        }
    }
    
    return(getPacksFromGRanges(collapsedRanges))
}

#' @title
#' packFinder Algorithm Pipeline
#'
#' @description
#' General use pipeline function for the Pack-TYPE transposon
#' finding algorithm.
#'
#' @param tirSeq
#' A \code{\link[Biostrings:DNAString-class]{DNAString}} 
#' object containing the TIR sequence to be searched for.
#'
#' @param Genome
#' A \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} 
#' object to be searched.
#'
#' @param mismatch
#' The maximum edit distance to be considered for TIR 
#' matches (indels + substitions). See 
#' \code{\link[Biostrings]{matchPattern}} for details.
#'
#' @param elementLength
#' The maximum element length to be considered, as a vector 
#' of two integers. E.g. \code{c(300, 3500)}
#'
#' @param tsdLength
#' Integer referring to the length of the flanking TSD region.
#' 
#' @param tsdMismatch
#' An integer referring to the allowable mismatch 
#' (substitutions or indels) between a transposon's TSD
#' sequences. \code{\link[Biostrings]{matchPattern}} from Biostrings 
#' is used for pattern matching. 
#'
#' @details
#' Finds potential pack-TYPE elements based on:
#' \itemize{
#'     \item Similarity of TIR sequence to \code{tirSeq}
#'     \item Proximity of potential TIR sequences
#'     \item Directionality of TIR sequences
#'     \item Similarity of TSD sequences
#' }
#'
#' The algorithm finds potential forward and reverse TIR 
#' sequences using \code{\link{identifyTirMatches}} and 
#' their associated TSD sequence via \code{\link{getTsds}}. 
#' The main filtering stage, 
#' \code{\link{identifyPotentialPackElements}}, filters 
#' matches to obtain a dataframe of potential PACK elements. 
#' Note that this pipeline does not consider the 
#' possibility of discovered elements being autonomous 
#' elements, so it is recommended to cluster and/or BLAST 
#' elements for further analysis. Furthermore, only exact TSD 
#' matches are considered, so supplying long sequences for 
#' TSD elements may lead to false-negative results.
#'
#' @return
#' A dataframe, containing elements 
#' identified by thealgorithm. These may be autonomous or 
#' pack-TYPE elements. Will contain the following features:
#' \itemize{
#'     \item start - the predicted element's start base 
#'     sequence position.
#'     \item end - the predicted element's end base 
#'     sequence position.
#'     \item seqnames - character string referring to the 
#'     sequence name in \code{Genome} to which \code{start} 
#'     and \code{end} refer to.
#'     \item width - the width of the predicted element.
#'     \item strand - the strand direction of the 
#'     transposable element. This will be set to "*" as the 
#'     \code{packSearch} function does not consider 
#'     transposons to have a direction - only TIR sequences. 
#'     Passing the \code{packMatches} dataframe to 
#'     \code{\link{packClust}} will assign a direction to 
#'     each predicted Pack-TYPE element.
#' }
#' 
#' This dataframe is in the format produced by 
#' coercing a \code{link[GenomicRanges:GRanges-class]{GRanges}} 
#' object to a dataframe: \code{data.frame(GRanges)}. Downstream 
#' functions, such as \code{\link{packClust}}, use this 
#' dataframe to manipulate predicted transposable elements.
#'
#' @note
#' This algorithm does not consider:
#' \itemize{
#'     \item Autonomous elements - autonomous elements will 
#'     be predicted by this algorithm as there is no BLAST 
#'     step. It is recommended that, after clustering 
#'     elements using \code{\link{packClust}}, the user 
#'     analyses each group to determine which predicted 
#'     elements are autonomous and which are likely 
#'     Pack-TYPE elements. Alternatively, databases such as 
#'     Repbase (\url{https://www.girinst.org/repbase/}) 
#'     supply annotations for autonomous transposable 
#'     elements that can be used to filter autonomous matches.
#'     \item TSD Mismatches - if two TIRs do not have exact 
#'     matches for their terminal site duplications they 
#'     will be ignored. Supplying longer TSD sequences will 
#'     likely lead to a lower false-positive rate, however 
#'     may also cause a greater rate of false-negative results.
#' }
#' 
#' Pattern matching is done via \code{\link[Biostrings]{matchPattern}}.
#'
#' @seealso
#' \code{\link{identifyTirMatches}}, \code{\link{getTsds}},
#' \code{\link{identifyPotentialPackElements}}, \code{\link{packClust}},
#' \code{\link{packMatches}}, 
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}},
#' \code{\link[Biostrings:DNAString-class]{DNAString}},
#' \code{\link[Biostrings]{matchPattern}}
#' 
#' 
#' @examples 
#' data(arabidopsisThalianaRefseq)
#'
#' packMatches <- packSearch(
#'     Biostrings::DNAString("CACTACAA"),
#'     arabidopsisThalianaRefseq,
#'     elementLength = c(300, 3500),
#'     tsdLength = 3
#' )
#' 
#' @author
#' Jack Gisby
#' 
#' @export

packSearch <- function(tirSeq, Genome, mismatch = 0, 
                        elementLength, tsdLength, tsdMismatch = 0) {
    
    searchCheck(mismatch, tsdLength, elementLength, tirSeq, Genome)

    message("Getting forward matches")
    forwardMatches <- identifyTirMatches(tirSeq = tirSeq, Genome = Genome,
        mismatch = mismatch, strand = "+", tsdLength = tsdLength)
    forwardMatches$TSD <- getTsds(tirMatches = forwardMatches, Genome = Genome,
        tsdLength = tsdLength, strand = "+")

    message(length(forwardMatches[, 1]), " forward matches identified.")
    message("Getting reverse matches")
    reverseMatches <- identifyTirMatches(
        tirSeq = Biostrings::reverseComplement(tirSeq), Genome = Genome, 
        mismatch = mismatch, strand = "-", tsdLength = tsdLength
    )
    reverseMatches$TSD <- getTsds(tirMatches = reverseMatches, Genome = Genome,
        tsdLength = tsdLength, strand = "-")

    message(length(reverseMatches[, 1]), " reverse matches identified.")
    if (length(forwardMatches[, 1]) == 0 | length(reverseMatches[, 1]) == 0) {
        message("No matches identified")
        return(NULL)
    }

    message("Filtering matches based on TSD sequences")
    packMatches <- identifyPotentialPackElements(
        forwardMatches = forwardMatches, reverseMatches = reverseMatches,
        Genome = Genome, elementLength = elementLength, 
        tsdMismatch = tsdMismatch
    )
    packMatches$TSD <- getTsds(tirMatches = packMatches, Genome = Genome,
                                tsdLength = tsdLength, strand = "+")
    
    message("Initial filtering complete. ", length(packMatches[, 1]), 
            " elements predicted.")
    return(packMatches)
}

searchCheck <- function(mismatch, tsdLength, elementLength, tirSeq, Genome) {
    if (!is.numeric(mismatch) | !is.numeric(tsdLength)) {
        stop("Arguments 'mismatch' and 'tsdLength' must be integers")
    }
    if (!is.vector(elementLength) | length(elementLength) != 2) {
        stop("Argument 'elementLength' must be a vector of minimum 
            and maximum transposon lengths")
    }
    if (!is.numeric(elementLength[1]) | !is.numeric(elementLength[2])) {
        stop("Vector 'elementLength' must contain integers")
    }
    if (!methods::is(tirSeq, "DNAString")) {
        stop("Argument 'tirSeq' must be of type 
            Biostrings::DNAString or character")
    }
    if (!methods::is(Genome, "DNAStringSet")) {
        stop("Argument 'Genome' must be of type Biostrings::DNAStringSet.
            You may convert files using Biostrings::readDNAStringSet
            or convert objects using Biostrings::DNAStringSet")
    }
}

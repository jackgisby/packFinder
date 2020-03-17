#' @title
#' Pipeline for BLAST/Classification of PackTYPE Elements
#'
#' @description
#' Run BLAST against user-specified databases of 
#' non-transposon and transposon-relates proteins.
#' Can be used to classify transposons based on 
#' their internal sequences. 
#'
#' @param packMatches
#' A dataframe of potential Pack-TYPE transposable elements, 
#' in the format given by \code{\link{packSearch}}. This 
#' dataframe is in the format produced by coercing a 
#' \code{link[GenomicRanges:GRanges-class]{GRanges}} 
#' object to a dataframe: \code{data.frame(GRanges)}. 
#' Will be saved as a FASTA file for VSEARCH.
#'
#' @param Genome
#' A DNAStringSet object containing sequences referred to 
#' in \code{packMatches} (the object originally used to 
#' predict the transposons \code{\link{packSearch}}). 
#' 
#' @param blastPath
#' Path to the BLAST+ executable, or name of 
#' the BLAST+ application for Linux/MacOS users.
#' 
#' @param protDb
#' For assigning Pack-TYPE elements. 
#' Path to the blast database containing nucleotide or protein
#' sequences to be matched against internal transposon 
#' sequences. Can be generated 
#' using BLAST+, or with 
#' \code{link{makeBlastDb}}.
#' 
#' @param autoDb
#' For assigning autonomous elements. 
#' Path to the blast database containing nucleotide or protein
#' sequences to be matched against internal transposon 
#' sequences. Can be generated 
#' using BLAST+, or with 
#' \code{link{makeBlastDb}}.
#' 
#' @param minE
#' Blast results with e values greater than
#' the specified cutoff will be ignored. This will 
#' be passed to BLASTN and applied to both transposon 
#' and non-transposon matches.
#' 
#' @param blastTask
#' Type of BLAST+ task, defaults to "blastn-short".
#' 
#' @param maxHits
#' Maximum hits returned by BLAST+ per query.
#' 
#' @param threads
#' Allowable number of threads to be utilised by BLAST+.
#' 
#' @param saveFolder
#' Directory to save BLAST+ results in; defaults 
#' to the working directory.
#' 
#' @param tirCutoff
#' How many bases to ignore at the terminal ends of the 
#' transposons to prevent hits to TIR sequences.
#' 
#' @param autoCutoff
#' Blast results for transposon-related elements will be 
#' filtered to 
#' ignore those with e values above the specified cutoff.
#' 
#' @param autoLength
#' Blast results for transposon-related elements containing 
#' hits with 
#' alignment lengths lower than this value will be ignored
#' 
#' @param autoIdentity
#' Blast results for transposon-related elements containing 
#' hits with 
#' sequence identities lower than this value will be ignored
#' 
#' @param autoScope
#' If specified, transposon-related blast results below the 
#' specified value
#' will be ignored. Note that the dataframe of transposon
#' matches must also be supplied to calculate scope. Scope is 
#' the proportion of the transposon's internal sequence 
#' occupied by the BLAST hit. 
#' 
#' @param protCutoff
#' Blast results for genic/other matches will be 
#' filtered to 
#' ignore those with e values above the specified cutoff.
#' 
#' @param protLength
#' Blast results for genic/other matches containing 
#' hits with 
#' alignment lengths lower than this value will be ignored
#' 
#' @param protIdentity
#' Blast results for genic/other matches containing 
#' hits with 
#' sequence identities lower than this value will be ignored
#' 
#' @param protScope
#' If specified, genic/other blast matches below the 
#' specified value
#' will be ignored. Note that the dataframe of transposon
#' matches must also be supplied to calculate scope. Scope is 
#' the proportion of the transposon's internal sequence 
#' occupied by the BLAST hit. 
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
#' \code{\link{blastAnalysis}}, \code{\link{packSearch}},
#' \code{\link{readBlast}}, \code{\link{blastAnnotate}}
#' 
#' @examples
#' \dontrun{
#' packMatches <- data(packMatches)
#' Genome <- data(arabidopsisThalianaRefseq)
#' 
#' packBlast(packMatches, Genome, 
#'     protDb = "C:/data/TAIR10_CDS", 
#'     autoDb = "C:/data/TAIR10_transposons", 
#'     blastPath = "C:/blast/bin/blastn.exe")
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
                          minE = autoCutoff, length = autoLength, 
                          identity = autoIdentity, removeExactMatches = FALSE, 
                          scope = autoScope, packMatches = packMatches)
    
    protHits <- readBlast(file.path(saveFolder, "protHits.blast"), 
                          minE = protCutoff, length = protLength, 
                          identity = protIdentity, removeExactMatches = TRUE,
                          scope = protScope, packMatches = packMatches)
    
    return(blastAnnotate(protHits, autoHits, packMatches))
}

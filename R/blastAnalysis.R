#' @title
#' BLAST Analysis of PackTYPE Elements 
#'
#' @description
#' Run BLAST against user-specified databases of 
#' non-transposon and transposon-related proteins.
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

blastAnalysis <- function(packMatches, Genome, blastPath, protDb = NULL, 
                          autoDb = NULL, minE = 1e-3, 
                          blastTask = "blastn-short", maxHits = 100, 
                          threads = 1, saveFolder = NULL, tirCutoff = 0) {
    
    if (is.null(saveFolder)) {
        saveFolder <- getwd()
    }
    
    protDb <- file.path(getwd(), protDb)
    autoDb <- file.path(getwd(), autoDb)
    
    packMatches$id <- rownames(packMatches)
    alternatePacksToFasta(packMatches, 
                          file.path(saveFolder, "packMatches.fasta"), 
                          Genome, tirCutoff)
    fastaFile <- file.path(saveFolder, "packMatches.fasta")
    
    if (!is.null(autoDb)) {
        executeBlast("auto", fastaFile, autoDb,
                     minE, maxHits, threads, blastPath, saveFolder, blastTask)
    }
    
    if (!is.null(protDb)) {
        executeBlast("prot", fastaFile, protDb,
                     minE, maxHits, threads, blastPath, saveFolder, blastTask)
    }
}

executeBlast <- function(type, fastaFile, blastDb, minE, maxHits, 
                         threads, blastPath, saveFolder, blastTask) {
    
    system2(blastPath, args = paste0(
        "-db ", blastDb, " -query ", fastaFile,
        " -out ", file.path(saveFolder, paste0(type, "Hits.blast")),
        " -evalue ", minE,
        " -max_target_seqs ", maxHits," -num_threads ", threads, 
        " -task ", blastTask, " -word_size 7", " -outfmt 6"
    ))
}

alternatePacksToFasta <- function(packMatches, file, Genome, i) {
    
    file.create(file)
    
    # create package-specific FASTA file
    for (match in seq_len(length(packMatches[, 1]))) {
        seq <- Genome[names(Genome) == packMatches$seqnames[match]][[1]]
        seq <- seq[(packMatches$start[match] + i):(packMatches$end[match] - i)]
        write(c(
            paste0(
                ">", packMatches$id[match],
                " | seq = ", packMatches$seqnames[match],
                " | start = ", packMatches$start[match],
                " | end = ", packMatches$end[match],
                " | width = ", packMatches$width[match],
                " | strand = ", packMatches$strand[match],
                " | TSD = ", packMatches$TSD[match]
            ),
            as.character(seq)
        ),
        file = file,
        append = TRUE
        )
    }
}


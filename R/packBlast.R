#' @title
#' BLAST Analysis of PackTYPE Elements 
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
#'     proteinDb = ", 
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

packBlast <- function(packMatches, Genome, proteinDb, autoDb, blastPath,
                      minE = 1e-3, blastTask = "blastn-short",
                      maxHits = 100, threads = 1, saveFolder = NULL,
                      tirCutoff = 0) {
    if (is.null(saveFolder)) {
        saveFolder <- getwd()
    }
    
    proteinDb <- file.path(getwd(), proteinDb)
    autoDb <- file.path(getwd(), autoDb)
    
    packMatches$id <- rownames(packMatches)
    alternatePacksToFasta(packMatches, file.path(saveFolder, "packMatches.fasta"), 
                          Genome, tirCutoff)
    fastaFile <- file.path(saveFolder, "packMatches.fasta")
    
    executeBlast("auto", fastaFile, autoDb, 
              minE, maxHits, threads, blastPath, saveFolder, blastTask)
    
    executeBlast("prot", fastaFile, proteinDb, 
                 minE, maxHits, threads, blastPath, saveFolder, blastTask)
}

executeBlast <- function(type, fastaFile, blastDb, minE, maxHits, 
                      threads, blastPath, saveFolder, blastTask) {
    
    system2(blastPath, args = paste0(
        "-db ", blastDb," -query ", fastaFile,
        " -out ", file.path(saveFolder, paste0(type, "Blast.blast")),
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


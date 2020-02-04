#' @title
#' Make Blast Database 
#'
#' @description
#' Generates a BLAST database to be queried. Required 
#' for identifying sequences using the BLAST+ 
#' software. 
#'
#' @param fastaFile
#' FASTA file containing sequences to generate a BLAST 
#' database from. 
#' 
#' @param dbPath
#' Path to save the BLAST database to. 
#' 
#' @param blastPath
#' Path/name of BLAST program to use. Name of
#' the application for Linux/MacOS, absolute 
#' path for the executable for windows users. 
#' 
#' @param dbType
#' Type of BLAST database to create, e.g. "nucl"
#' for a nucleotide database. 
#'     
#' @seealso 
#' \code{\link{packSearch}}
#' 
#' @examples
#' \dontrun{
#' makeBlastDb("genes.fasta", "blastdb.db", "C:/blast.exe")
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

makeBlastDb <- function(fastaFile, dbPath, blastPath, dbType = "nucl") {
    system2(blastPath, args = paste0(
        "-in ", file.path(getwd(), fastaFile), " ", 
        "-input_type fasta -dbtype ", dbType, 
        " -parse_seqids -out ", dbPath, " ",
        "-title packfinder_DB"   
    ))
}

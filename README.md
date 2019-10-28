# packFinder
A package for the de novo discovery of pack-TYPE transposons. Transposons are detected in DNA sequences based on conserved terminal inverted repeat sequences and presence of terminal site duplications. packFinder allows users to search a given Genome for Pack-TYPE transposons with minimal input and setup.

![**Important structural features of Pack-TYPE transposons**](vignettes/tirSeq.jpg)

## Bioconductor Submission
R CMD check was run in <2 minutes with no warnings, errors or notes. BiocCheck was run with a single note: 
"Usage of dontrun{} / donttest{} found in man page examples."

Two functions, packClust and packAlign, are not run in examples or vignettes. These functions use the command line tool VSEARCH (https://github.com/torognes/vsearch) to cluster putative transposable elements - could not find an equivalent tool in R that clustered sequences with the same performance and accuracy. Note this is not part of the main package annotation pipeline; the main functions, such as packSearch, can be run without installation of VSEARCH and are run in examples/vignettes in the package. 

## Using packFinder
The packFinder vignette includes a full walkthrough of the package. To get started quickly, it is easiest to download a genome of interest using the biomartr package and follow the steps outlined in the vignette. Additionally, VSEARCH installation instructions can be found here.

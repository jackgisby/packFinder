# packFinder <img src="inst/packFinder_hex.png" align="right" height="174" width="150" />

A package for the de novo discovery of pack-TYPE transposons. Transposons are detected in DNA sequences based on conserved terminal inverted repeat sequences and presence of terminal site duplications. packFinder allows users to search a given Genome for Pack-TYPE transposons with minimal input and setup.

## Installation
packFinder is available from the development branch of Bioconductor (https://bioconductor.org/packages/packFinder/). With R 4.0.0 installed, the following code can be used to install the package:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("packFinder")
```

Additionally, in order to use the functions packClust and packAlign, the VSEARCH command line tool must be installed. For Linux and MacOS systems, correct installation of VSEARCH should allow users to use all functions within packFinder; for windows users, the absolute path to the VSEARCH executable file must be specified when calling packFinder clustering and alignment functions. See the vignette, documentation or the VSEARCH github for further information (https://github.com/torognes/vsearch).

## Using packFinder
The packFinder vignette includes a full walkthrough of the package. To get started quickly, it is easiest to download a genome of interest using the biomartr package and follow the steps outlined in the vignette ().

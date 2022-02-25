# packFinder <img src="inst/packFinder_hex.png" align="right" height="174" width="150" />

A package for the *de novo* discovery of Pack-TYPE transposons. Transposons are detected in DNA sequences based on conserved terminal inverted repeat sequences and presence of terminal site duplications. `packFinder` allows users to search a given Genome for Pack-TYPE transposons with minimal input and setup. For a more in-depth explanation, see our [paper in PLOS Genetics](https://doi.org/10.1371/journal.pgen.1010078).

## Installation
`packFinder` is available from the [Bioconductor project](https://bioconductor.org/packages/packFinder/). With R >=4.0.0 installed, the following code can be used to install the package:
```
if (!require("BiocManager"))
    install.packages("BiocManager")
    
BiocManager::install("packFinder")
```

In order to use the functions `packClust` and `packAlign`, the VSEARCH command line tool must be installed. For Linux and MacOS systems, correct installation of VSEARCH should allow users to use all functions within `packFinder`; for windows users, the absolute path to the VSEARCH executable file must be specified when calling `packFinder` clustering and alignment functions. See the vignette, documentation or the [VSEARCH github](https://github.com/torognes/vsearch) for further information.

Additionally, utilites are provided with the package for the functional annotation of internal transposon sequences using BLAST+; this application can be installed from [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

## Using packFinder
The `packFinder` vignette includes a full walkthrough of the package. To get started quickly, it is easiest to download a genome of interest using the `biomartr` package and follow the steps outlined in [the vignette](https://bioconductor.org/packages/release/bioc/vignettes/packFinder/inst/doc/packFinder.html).

## Citation
If you use `packFinder` in your work, please cite the following publication:

> Gisby, J. S. & Catoni, M. (2022) [The widespread nature of Pack-TYPE transposons reveals their importance for plant genome evolution.](https://doi.org/10.1371/journal.pgen.1010078) *PLOS Genetics* 18(2): e1010078.
> 

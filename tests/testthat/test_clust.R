context("clustering and alignment")

dir.create("tempTestOutput")
data("packMatches")
data("arabidopsisThalianaRefseq")

# cluster testing will be skipped if location of vsearch is not specified
# vSearchLocation <- "E:/vsearch-2.14.1-win-x86_64/vsearch.exe"
vSearchLocation <- NULL

if (!is.null(vSearchLocation)) {
    packClusts <- packClust(packMatches, arabidopsisThalianaRefseq, 
                            saveFolder = "data-raw/output", 
                            vSearchPath = vSearchLocation)
    
    packAlign <- packAlign(packMatches, arabidopsisThalianaRefseq, 
                           saveFolder = "data-raw/output", 
                           vSearchPath = vSearchLocation)
}

test_that("Clusters identified are as expected", {
    skip_if(is.null(vSearchLocation), 
            "Cluster functions cannot be tested without VSEARCH.")
    expect_equal(packClusts, packMatches)
    expect_equal(packAlign, 
                readUc("data-raw/output/vSearchPairwiseAlignment.uc", 
                        output = "alignment"))
})

unlink("tempTestOutput", recursive = TRUE)

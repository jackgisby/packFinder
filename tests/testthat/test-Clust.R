context("clustering and alignment")

data("packMatches")
data("arabidopsisThalianaRefseq")

if(!is.null(vSearchLocation)) {
  packClusts <- packClust(packMatches, arabidopsisThalianaRefseq, saveFolder = "data-raw/output/", vSearchPath = vSearchLocation)
  packAlign <- packAlign(packMatches, arabidopsisThalianaRefseq, saveFolder = "data-raw/output/", vSearchPath = vSearchLocation)
}

test_that("Clusters identified are as expected", {
  expect_equal(packClusts, packMatches)
  expect_equal(packAlign, readUc("data-raw/output/vSearchPairwiseAlignment.uc", output = "alignment"))
})

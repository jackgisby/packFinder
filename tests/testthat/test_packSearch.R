context("packSearch")

data("arabidopsisThalianaRefseq")
load("data-r/packClusts.rda")
subSeq <- Biostrings::DNAString("CACTACAA")

packMatches <- packSearch(subSeq,
  arabidopsisThalianaRefseq,
  mismatch = 0,
  elementLength = c(300, 3500),
  tsdLength = 3
)

test_that("dimensions of packMatches are correct", {
  expect_equal(nrow(packMatches), 29)
  expect_equal(ncol(packMatches), 6)

  expect_equal(nrow(packClusts), 29)
  expect_equal(ncol(packClusts), 7)
})

test_that("packSearch returns correct results", {
  expect_equal(packMatches$seqnames, packClusts$seqnames)
  expect_equal(packMatches$start, packClusts$start)
  expect_equal(packMatches$end, packClusts$end)
  expect_equal(packMatches$width, packClusts$width)
  expect_equal(packMatches$seq, packClusts$seq)
  expect_equal(packMatches$TSD, packClusts$TSD)
})
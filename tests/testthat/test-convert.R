context("convert")

data("packMatches")

packsToCsv(packMatches, file = "data-raw/output/packMatches.csv")
packsToFasta(packMatches, file = "data-raw/output/packMatches.fasta", Genome = arabidopsisThalianaRefseq)
packsGRanges <- packsToGRanges(packMatches)

getPacksFromCsv <- getPacksFromCsv(file = "data-raw/output/packMatches.csv")
getPacksFromFasta <- getPacksFromFasta("data-raw/output/packMatches.fasta")
getPacksFromGRanges <- getPacksFromGRanges(packsGRanges)

test_that("Conversion functions create correct output from sample data", {
  expect_equal(packsFromCsv, packMatches)
  expect_equal(packsFromFasta, packMatches)
  expect_equal(packsFromGRanges, packMatches)
})

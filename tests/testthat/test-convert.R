context("convert")

data("packMatches")

packsToCsv(packMatches, "data-raw/output/packMatches.csv")
packsToFasta(packMatches, "data-raw/output/packMatches.fasta")
packsGRanges <- packsToGRanges(packMatches)

packsFromCsv <- getPacksFromCsv("data-raw/output/packMatches.csv")
packsFromFasta <- getPacksFromFasta("data-raw/output/packMatches.fasta")
packsFromGRanges <- getPacksFromGRanges(packsGRanges)

test_that("Conversion functions create correct output from sample data", {
  expect_equal(packsFromCsv, packMatches)
  expect_equal(packsFromFasta, packMatches)
  expect_equal(packsFromGRanges, packMatches)
})

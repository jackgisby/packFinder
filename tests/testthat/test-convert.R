context("convert")

data("packMatches")

packsToCsv(packMatches, file = "data-raw/output/packMatches.csv")
packsToFasta(packMatches, file = "data-raw/output/packMatches.fasta", Genome = arabidopsisThalianaRefseq)
packsGRanges <- packsToGRanges(packMatches)

packsFromCsv <- getPacksFromCsv(file = "data-raw/output/packMatches.csv")
packsFromFasta <- getPacksFromFasta("data-raw/output/packMatches.fasta")
packsFromGRanges <- getPacksFromGRanges(packsGRanges)

#ignore differences in column type
packsFromCsv$seqnames <- as.factor(packsFromCsv$seqnames)
packsFromFasta$seqnames <- as.factor(packsFromFasta$seqnames)
packsFromFasta$start <- as.integer(packsFromFasta$start)
packsFromFasta$end <- as.integer(packsFromFasta$end)
packsFromFasta$width <- as.integer(packsFromFasta$width)
packsFromGRanges$strand <- as.character(packsFromGRanges$strand)

test_that("Conversion functions create correct output from sample data", {
  expect_equal(as.character(packsFromCsv), as.character(packMatches))
  expect_equal(packsFromFasta, subset(packMatches, select = -c(cluster)))
  expect_equal(packsFromGRanges, packMatches)
})

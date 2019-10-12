context("Conversion functions")

ucClust <- readUc("data-raw/clustOutput.uc")
ucAlign <- readUc("data-raw/alignOutput.uc")
blast6Out <- readBlast6Out("data-raw/clustOutput.blast6out")

test_that("object type is as expected", {
  expect_type(ucClust, "list")
  expect_type(ucAlign, "list")
  expect_type(blast6Out, "list")
})

test_that("object dimensions are as expected", {
  expect_equal(nrow(ucClust), 11)
  expect_equal(ncol(ucClust), 8)

  expect_equal(nrow(ucAlign), 16)
  expect_equal(ncol(ucAlign), 8)

  expect_equal(nrow(blast6Out), 1)
  expect_equal(ncol(blast6Out), 12)
})

test_that("object output is as expected", {
  expect_equal(max(ucClust$cluster), 4)
})

unlink("data-raw/output/*")

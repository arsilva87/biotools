library(testthat)
context("test tocher")

test_that("tocher works", {
  to <- tocher(dist(rnorm(10)))
  expect_output(str(to), "List of 7")
})

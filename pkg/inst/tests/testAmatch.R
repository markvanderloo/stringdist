
context("amatch: osa")

test_that("simple test and multiple edge cases",{
  expect_equal(amatch("aa",c("ba","bb"), method="osa"), 1L)
  expect_equal(amatch("","", method="osa"), 1L)
  expect_equal(amatch(NA,"a", method="osa"), NA_integer_)
  expect_equal(amatch(NA,"a", method="osa", nomatch=0L),0L)
  expect_equal(amatch(NA,NA, method="osa"), 1L)
  expect_equal(amatch(NA,NA, method="osa",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="osa",matchNA=FALSE, nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="osa",matchNA=FALSE, nomatch=7L), 7L)
  expect_equal(amatch("aa","bb", method="osa",maxDist=1), NA_integer_)
  expect_equal(amatch("aa","bb", method="osa",maxDist=1), NA_integer_)
})


context("amatch: Damerau-Levenshtein")

test_that("simple test and multiple edge cases",{
  expect_equal(amatch("aa", c("ba","bb"), method="dl"), 1L)
  expect_equal(amatch("","", method="dl"), 1L)
  expect_equal(amatch(NA,"a", method="dl"), NA_integer_)
  expect_equal(amatch(NA,"a", method="dl",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="dl"), 1L)
  expect_equal(amatch(NA,NA, method="dl",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="dl",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="dl",matchNA=FALSE,nomatch=7L), 7L)
  expect_equal(amatch("aa","bb", method="dl",maxDist=1), NA_integer_)
  expect_equal(amatch("aa","bb", method="dl",maxDist=1), NA_integer_)
})

context("amatch: Hamming")

test_that("simple test and multiple edge cases",{
  expect_equal(amatch("aa", c("ba","bb"), method="hamming"), 1L)
  expect_equal(amatch(NA,c(NA,NA),method="hamming"),1L)
  expect_equal(amatch("","", method="hamming"), 1L)
  expect_equal(amatch(NA,"a", method="hamming"), NA_integer_)
  expect_equal(amatch(NA,"a", method="hamming",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="hamming"), 1L)
  expect_equal(amatch(NA,NA, method="hamming",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="hamming",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="hamming",matchNA=FALSE,nomatch=7L), 7L)
  expect_equal(amatch("aa","bb", method="hamming",maxDist=1), NA_integer_)
})


#context("amatch: Jaro and Jaro-Winkler")
#
#test_that("simple test and multiple edge cases",{
#  expect_equal(amatch("aa", c("ba","bb"), method="jw"), 1L)
#  expect_equal(amatch(NA,c(NA,NA),method="jw"),1L)
#  expect_equal(amatch("","", method="jw"), 1L)
#  expect_equal(amatch(NA,"a", method="jw"), NA_integer_)
#  expect_equal(amatch(NA,"a", method="jw",nomatch=0L), 0L)
#  expect_equal(amatch(NA,NA, method="jw"), 1L)
#  expect_equal(amatch(NA,NA, method="jw",matchNA=FALSE), NA_integer_)
#  expect_equal(amatch(NA,NA, method="jw",matchNA=FALSE,nomatch=0L), 0L)
#  expect_equal(amatch(NA,NA, method="jw",matchNA=FALSE,nomatch=7L), 7L)
#  expect_equal(amatch("aa","bb", method="jw",maxDist=1), NA_integer_)
#})




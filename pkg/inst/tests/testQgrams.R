
context("qgrams")

test_that("qgram edge cases",{
  expect_equivalent(qgrams('a' ,      q=1), as.matrix(c(a=1)))   # basic test
  expect_equivalent(qgrams('aa',      q=1), as.matrix(c(a=2)))   # idem
  expect_equivalent(qgrams(c('a','a'),q=1), as.matrix(c(a=2)))   # count unique n-grams
  expect_equivalent(qgrams(c(NA,'a'), q=1), as.matrix(c(a=1)))   # skip NA's
  expect_equivalent(qgrams(NA,q=1), matrix(0,nrow=1,ncol=0))          # skip all
  expect_equivalent(qgrams(c("a","ab"), q=2), as.matrix(table("ab"))) # skip q>nchar
  expect_equivalent(qgrams(c("a"),q=2), matrix(0,nrow=1,ncol=0))      # skip all
  expect_equivalent(qgrams(c(''),q=0), as.matrix(table('')))          # empty string, q=0
})



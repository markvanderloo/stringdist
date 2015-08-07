
context("seq_dist")

# A simple test to see that everything is passed on to the correct
# algorithm
test_that("Methods are selected and computed correctly", {
  expect_equal(
    seq_dist(a = list(c(1L,2L,3L)), b = list(c(2L,1L,3L)), method="osa")
    , 1 )
  expect_equal(
    seq_dist(a = list(c(1L,2L,3L)), b = list(c(2L,1L,3L)), method="lv")
    , 2 )
  # the case setting 'dl' apart from 'osa'
  expect_equal(
    seq_dist(a = list(c(2L,1L)), b = list(c(1L,3L,2L)), method="dl")
    , 2 )
  expect_equal(
    seq_dist(a = list(c(1L,2L,3L)), b = list(c(1L,0L,3L)), method="hamming")
    , 1 )
  expect_equal(
    seq_dist(a = list(c(1L,2L,3L)), b = list(c(1L,0L,3L)), method="lcs")
    , 2 )
  expect_equal(
    seq_dist(a = list(c(1L,2L,3L)), b = list(c(1L,0L,3L)), method="qgram",q=2)
  , 4 )
  
  expect_equal(
    round(1-seq_dist(list(utf8ToInt("martha")),list(utf8ToInt("marhta")),method='jw'),3)
    , 0.944
  )
  expect_error(
    seq_dist(a = list(c(1L,2L,3L)), b = list(c(1L,0L,3L)), method="soundex")
  )
})  

test_that("Some edge cases",{
  expect_equal(length(seq_dist(list(),list(c(1L)))),0)
  expect_equal(length(seq_dist(list(),list())),0)
})

test_that("Elementary tests on seq_distmatrix",{
  expect_error(seq_distmatrix(1:10))
  expect_error(seq_distmatrix(1:10,list(1:10)))
  expect_error(seq_distmatrix(list(1:10),1:10))
  expect_equivalent(
    as.matrix(seq_distmatrix(list(1:3,2:4)) )
    , matrix(c(0,2,2,0),nrow=2)
  )
  expect_equal(
    as.matrix(seq_distmatrix(list(x=1:3,y=2:4),useNames="names") )
    , matrix(c(0,2,2,0),nrow=2,dimnames=list(c('x','y'),c('x','y')))
  )
  expect_equal(
    seq_distmatrix(list(x=1:3,y=2:4),list(x=1:3,y=2:4),useNames="names")
    , matrix(c(0,2,2,0),nrow=2,dimnames=list(c('x','y'),c('x','y')))
  )
  expect_equal(class(seq_distmatrix(list(1:3,2:4))),"dist")
  expect_equivalent(
    as.matrix(seq_distmatrix(list(1:3,2:4)),seq_distmatrix(list(1:3,2:4),list(1:3,2:4)) )
    , matrix(c(0,2,2,0),nrow=2)
  )
  
  
  
})



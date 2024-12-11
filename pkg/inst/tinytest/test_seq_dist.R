options(sd_num_thread=2)
## seq_dist
# tests against cases that used to segfault when we did not check
# NULL cases.
expect_error(seq_dist(a=list(c(1L,2L,3L)), b=list(c(2L,1L,3L)),nthread=1:4))
expect_error(seq_dist(a=list(c(1L,2L,3L)), b=list(c(2L,1L,3L)),nthread="foo"))
expect_error(seq_dist(a=list(c(1L,2L,3L)), b=list(c(2L,1L,3L)),nthread=integer(0)))
expect_error(seq_dist(a=list(c(1L,2L,3L)), b=list(c(2L,1L,3L)),nthread=NULL))

# A simple test to see that everything is passed on to the correct
# algorithm
## Methods are selected and computed correctly
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


## Conversion for non-integer-list arguments
  expect_equal(seq_dist(list(c(1,2,3)),list(c(2,3,4))),seq_dist(as.numeric(c(1,2,3)),as.numeric(c(2,3,4))))
  expect_equal(seq_dist(list(c(1,2,3)),list(c(2,3,4))),seq_dist(c(1,2,3), c(2,3,4)))
  expect_equal(seq_distmatrix(list(c(1,2,3)),list(c(2,3,4))), seq_distmatrix(as.numeric(c(1,2,3)),as.numeric(c(2,3,4))))  
  expect_equal(seq_distmatrix(list(c(1,2,3)),list(c(2,3,4))), seq_distmatrix(c(1,2,3),c(2,3,4)))  
  expect_equal(seq_distmatrix(list(c(1,2,3))),seq_distmatrix(c(1,2,3)))
  expect_equal(seq_distmatrix(list(c(1,2,3))),seq_distmatrix(as.numeric(c(1,2,3))))


## Some edge cases
  expect_equal(length(seq_dist(list(),list(c(1L)))),0)
  expect_equal(length(seq_dist(list(),list())),0)


## Elementary tests on seq_distmatrix

  expect_equivalent(seq_distmatrix(1:10),dist(0))
  expect_equivalent(seq_distmatrix(1:10,list(1:10)),matrix(0))
  expect_equivalent(
    as.matrix(seq_distmatrix(list(c(1,2,3),c(2,3,4))) )
    , matrix(c(0,2,2,0),nrow=2)
  )
  expect_equal(
    as.matrix(seq_distmatrix(list(x=c(1,2,3),y=c(2,3,4)),useNames="names") )
    , matrix(c(0,2,2,0),nrow=2,dimnames=list(c('x','y'),c('x','y')))
  )
  expect_equal(
    seq_distmatrix(list(x=c(1,2,3),y=c(2,3,4)),list(x=c(1,2,3),y=c(2,3,4)),useNames="names")
    , matrix(c(0,2,2,0),nrow=2,dimnames=list(c('x','y'),c('x','y')))
  )
  expect_equal(class(seq_distmatrix(list(c(1,2,3),c(2,3,4)))),"dist")
  expect_equivalent(
    as.matrix(seq_distmatrix(list(c(1,2,3),c(2,3,4))),seq_distmatrix(list(c(1,2,3),c(2,3,4)),list(c(1,2,3),c(2,3,4))) )
    , matrix(c(0,2,2,0),nrow=2)
  )
  



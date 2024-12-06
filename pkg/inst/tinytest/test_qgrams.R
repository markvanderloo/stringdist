options(sd_num_thread=2)
## qgrams

## qgram edge cases
  expect_equivalent(qgrams('a' ,      q=1), as.matrix(c(a=1)))   # basic test
  expect_equivalent(qgrams('aa',      q=1), as.matrix(c(a=2)))   # idem
  expect_equivalent(qgrams(c('a','a'),q=1), as.matrix(c(a=2)))   # count unique n-grams
  expect_equivalent(qgrams(c(NA,'a'), q=1), as.matrix(c(a=1)))   # skip NA's
  expect_equivalent(qgrams(NA,q=1), matrix(0,nrow=1,ncol=0))          # skip all
  expect_equivalent(qgrams(c("a","ab"), q=2), as.matrix(table("ab"))) # skip q>nchar
  expect_equivalent(qgrams(c("a"),q=2), matrix(0,nrow=1,ncol=0))      # skip all
  expect_equivalent(qgrams(c(''),q=0), matrix(table('')))          # empty string, q=0


## qgrams
  expect_equivalent(qgrams("a",q=1),array(1,dim=c(1,1)))
  expect_equivalent(qgrams("a",q=1,useBytes=TRUE),array(1,dim=c(1,1)))


## seq_qgrams
  expect_equivalent(
    seq_qgrams(1:3,2:4,q=2)
    ,matrix(c(
       1,2,1,0
      ,2,3,1,1
      ,3,4,0,1
    ),nrow=3,byrow=TRUE)
  )


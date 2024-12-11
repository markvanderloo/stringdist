options(sd_num_thread=2)
## amatch: Optimal String Alignment
# tests against cases that used to segfault when we did not check
# NULL cases.
expect_error(amatch("a","b",nthread=1:4))
expect_error(amatch("a","b",nthread="foo"))
expect_error(amatch("a","b",nthread=integer(0)))
expect_error(amatch("a","b",nthread=NULL))


## simple test and multiple edge cases
  expect_equal(amatch("aa",c("ba","bb"), method="osa",maxDist=1L), 1L)
  expect_equal(amatch("aa",c("bb","bb"), method="osa",maxDist=1L), NA_integer_)
  expect_equal(amatch("aa",c("bbb"), method="osa",maxDist=2L), NA_integer_)
  expect_equal(amatch("bbb",c("aa"), method="osa",maxDist=2L), NA_integer_)
  expect_equal(amatch("","", method="osa"), 1L)
  expect_equal(amatch(NA,"a", method="osa"), NA_integer_)
  expect_equal(amatch(NA,"a", method="osa", nomatch=0L),0L)
  expect_equal(amatch(NA,NA, method="osa"), 1L)
  expect_equal(amatch(NA,NA, method="osa",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="osa",matchNA=FALSE, nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="osa",matchNA=FALSE, nomatch=7L), 7L)
  expect_equal(amatch("aa","bb", method="osa",maxDist=1), NA_integer_)
  expect_equal(amatch("aa","bb", method="osa",maxDist=1), NA_integer_)
  expect_equal(amatch(c("m","fem"),c("male","female"),method="osa",maxDist=Inf), c(1,2))




## amatch: Damerau-Levenshtein

## simple test and multiple edge cases
  expect_equal(amatch("aa", c("ba","bb"), method="dl",maxDist=1L), 1L)
  expect_equal(amatch("aa",c("bb","bb"), method="dl",maxDist=1L), NA_integer_)
  expect_equal(amatch("aa",c("bbb"), method="dl",maxDist=2L), NA_integer_)
  expect_equal(amatch("bbb",c("aa"), method="dl",maxDist=2L), NA_integer_)
  expect_equal(amatch("","", method="dl"), 1L)
  expect_equal(amatch(NA,"a", method="dl"), NA_integer_)
  expect_equal(amatch(NA,"a", method="dl",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="dl"), 1L)
  expect_equal(amatch(NA,NA, method="dl",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="dl",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="dl",matchNA=FALSE,nomatch=7L), 7L)
  expect_equal(amatch("aa","bb", method="dl",maxDist=1), NA_integer_)
  expect_equal(amatch("aa","bb", method="dl",maxDist=1), NA_integer_)
  expect_equal(amatch(c("m","fem"),c("male","female"),method="dl",maxDist=Inf), c(1,2))


## amatch: Hamming

## simple test and multiple edge cases
  expect_equal(amatch("aa", c("ba","bb"), method="hamming",maxDist=1L), 1L)
  expect_equal(amatch("aa",c("bb","bb"), method="hamming",maxDist=1L), NA_integer_)
  expect_equal(amatch(NA,c(NA,NA),method="hamming"),1L)
  expect_equal(amatch("","", method="hamming"), 1L)
  expect_equal(amatch(NA,"a", method="hamming"), NA_integer_)
  expect_equal(amatch(NA,"a", method="hamming",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="hamming"), 1L)
  expect_equal(amatch(NA,NA, method="hamming",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="hamming",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="hamming",matchNA=FALSE,nomatch=7L), 7L)
  expect_equal(amatch("aa","bb", method="hamming",maxDist=1), NA_integer_)



## amatch: Jaro and Jaro-Winkler

## simple test and multiple edge cases
  expect_equal(amatch("aa", c("ba","bb"), method="jw",maxDist=1L), 1L)
  expect_equal(amatch("aa",c("bb","bb"), method="jw",maxDist=0.5), NA_integer_)
  expect_equal(amatch(NA,c(NA,NA),method="jw"),1L)
  expect_equal(amatch("","", method="jw"), 1L)
  expect_equal(amatch(NA,"a", method="jw"), NA_integer_)
  expect_equal(amatch(NA,"a", method="jw",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="jw"), 1L)
  expect_equal(amatch(NA,NA, method="jw",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="jw",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="jw",matchNA=FALSE,nomatch=7L), 7L)
  expect_equal(amatch(c("m","fem"),c("male","female"),method="jw",maxDist=Inf), c(1,2))

## amatch: Longest Common Substring

## simple test and multiple edge cases
  expect_equal(amatch("aa", c("ba","bb"), method="lcs",maxDist=2L), 1L)
  expect_equal(amatch("aa",c("bb","bb"), method="lcs",maxDist=1L), NA_integer_)
  expect_equal(amatch("aa",c("bbb"), method="lcs",maxDist=2L), NA_integer_)
  expect_equal(amatch("bbb",c("aa"), method="lcs",maxDist=2L), NA_integer_)
  expect_equal(amatch(NA,c(NA,NA),method="lcs"),1L)
  expect_equal(amatch("","", method="lcs"), 1L)
  expect_equal(amatch(NA,"a", method="lcs"), NA_integer_)
  expect_equal(amatch(NA,"a", method="lcs",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="lcs"), 1L)
  expect_equal(amatch(NA,NA, method="lcs",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="lcs",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="lcs",matchNA=FALSE,nomatch=7L), 7L)


## amatch: Levenshtein

## simple test and multiple edge cases
  expect_equal(amatch("aa", c("ba","bb"), method="lv",maxDist=1L), 1L)
  expect_equal(amatch("aa",c("bb","bb"), method="lv",maxDist=1L), NA_integer_)
  expect_equal(amatch("aa",c("bbb"), method="lv",maxDist=2L), NA_integer_)
  expect_equal(amatch("bbb",c("aa"), method="lv",maxDist=2L), NA_integer_)
  expect_equal(amatch(NA,c(NA,NA),method="lv"),1L)
  expect_equal(amatch("","", method="lv"), 1L)
  expect_equal(amatch(NA,"a", method="lv"), NA_integer_)
  expect_equal(amatch(NA,"a", method="lv",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="lv"), 1L)
  expect_equal(amatch(NA,NA, method="lv",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="lv",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="lv",matchNA=FALSE,nomatch=7L), 7L)
  expect_equal(amatch(c("m","fem"),c("male","female"),method="lv",maxDist=Inf), c(1,2))


## amatch: qgrams

## simple test and multiple edge cases
  expect_equal(amatch("aa", c("ba","bb"), method="qgram",maxDist=2), 1L)
  expect_equal(amatch("aa",c("bb","bb"), method="qgram",maxDist=1L), NA_integer_)
  expect_equal(amatch(NA,c(NA,NA),method="qgram"),1L)
  expect_equal(amatch("","", method="qgram", q=0), 1L)
  expect_equal(amatch(NA,"a", method="qgram"), NA_integer_)
  expect_equal(amatch(NA,"a", method="qgram",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="qgram"), 1L)
  expect_equal(amatch(NA,NA, method="qgram",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="qgram",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="qgram",matchNA=FALSE,nomatch=7L), 7L)
  # bug reported by Robert Carlson
  expect_equal( 
      amatch("600 EXAMPLE AVE NJ 8629", 
        c("2100 EXAMPLE AVE NJ 8619", "600 EXAMPLE AVE NJ 8629"), method="jaccard")
      , 2L)



## amatch: Soundex

## simple test and multiple edge cases
  expect_equal(amatch("smith", c("smyth","smelt"), method="soundex"), 1L)
  expect_equal(amatch("smith",c("bb","bb"), method="soundex"), NA_integer_)
  expect_equal(amatch("smith",c("whashington"), method="soundex"), NA_integer_)
  expect_equal(amatch("washington",c("smith"), method="soundex"), NA_integer_)
  expect_equal(amatch(NA,c(NA,NA),method="soundex"),1L)
  expect_equal(amatch("","", method="soundex"), 1L)
  expect_equal(amatch(NA,"a", method="soundex"), NA_integer_)
  expect_equal(amatch(NA,"a", method="soundex",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="soundex"), 1L)
  expect_equal(amatch(NA,NA, method="soundex",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="soundex",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="soundex",matchNA=FALSE,nomatch=7L), 7L)


## amatch: useBytes

## bytewise matching differs from character wise matching
  x <- paste0('Mot',intToUtf8(0x00F6),'rhead') 
  y <- c('bastard','Motorhead') 
  jwdist <- round(1-(1/3)*(8/9 + 8/10 + 1),3)

  expect_equal(amatch(x, y, method='dl', maxDist=2, useBytes=TRUE), 2);
  expect_equal(amatch(x, y, method='dl', maxDist=1, nomatch=0L, useBytes=TRUE), 0L);
  expect_equal(amatch(x, x, method='hamming',maxDist=1L, useBytes=TRUE),1); 
  expect_equal(amatch(x, y, method='hamming',nomatch=0L,useBytes=TRUE), 0L);
  expect_equal(amatch(x, y, method='jw', maxDist=1.0, useBytes=TRUE), 2);
  expect_equal(amatch(x, y, method='jw', maxDist=jwdist-0.01, nomatch=0L, useBytes=TRUE), 0L);
  expect_equal(amatch(x, y, method='lcs',maxDist=3, useBytes=TRUE),2L); 
  expect_equal(amatch(x, y, method='lcs',maxDist=2, useBytes=TRUE, nomatch=0L), 0L);
  expect_equal(amatch(x, y, method='lv',maxDist=2, useBytes=TRUE),2L); 
  expect_equal(amatch(x, y, method='lv',maxDist=1, useBytes=TRUE, nomatch=0L), 0L);
  expect_equal(amatch(x, y, method='osa',maxDist=2, useBytes=TRUE),2L); 
  expect_equal(amatch(x, y, method='osa',maxDist=1, useBytes=TRUE, nomatch=0L), 0L);
  expect_equal(amatch(x, y, method='qgram',maxDist=7, q=3, useBytes=TRUE),2L); 
  expect_equal(amatch(x, y, method='qgram',maxDist=6, q=3, useBytes=TRUE, nomatch=0L), 0L);


## seq_amatch

## Input checks for seq_amatch
  expect_equal(seq_amatch(list(1:10),list(1:10)),seq_amatch(1:10,1:10))
  expect_equal(seq_amatch(list(1:10),list(1:10)),seq_amatch(as.numeric(1:10),as.numeric(1:10)))


## Some elementary tests for seq_amatch and seq_ain
  
  x <- list(c(1L,3L,2L))
  table <- list(
    c(2L,3L,4L)
    ,c(1L,2L,3L)
  )
  expect_equal(seq_amatch(x,table,maxDist=3), 2)
  expect_equal(seq_amatch(list(NA_integer_),table,maxDist=3),NA_integer_ )
  expect_true(seq_ain(x,table,maxDist=3))
  expect_false(seq_ain(x,table))


expect_error(seq_amatch(x,table,nthread=1:4))
expect_error(seq_amatch(x,table,nthread="foo"))
expect_error(seq_amatch(x,table,nthread=integer(0)))
expect_error(seq_amatch(x,table,nthread=NULL))



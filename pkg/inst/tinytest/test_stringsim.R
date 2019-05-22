options(sd_num_thread=2)
## stringsim

# We expect that two completely different strings have a similarity of
# 0 and two completely equal strings a similarity of 1
methods <- c("osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex")
for (method in methods) {

    expect_equal(stringsim("bb", "cc", method=method), 0)
    expect_equal(stringsim("bb", "bb", method=method), 1)
}


## edgecases
for (method in methods[c(1:5,9:10)]){
    expect_equal(stringsim(c("a", ""), "", method=method), c(0, 1))
    
    expect_equal(stringsim(c("kkk", "bbb"), "bbb", method=method), 
      stringsim("bbb", c("kkk", "bbb"), method=method))
}
for (method in methods[6:8]){
    expect_equal(stringsim(c("a", ""), "", method=method,q=0), c(1, 1))
    
    expect_equal(stringsim(c("kkk", "bbb"), "bbb", method=method,q=0), 
                stringsim("bbb", c("kkk", "bbb"), method=method,q=0))
}




## Stringsim gets correct values with or without useBytes
  x <- "ao"
  y <- paste0("a",intToUtf8(0x00F6)) # o-umlaut
  expect_equal(stringsim(x,y,method="osa", useBytes=FALSE), 1-1/2)
  expect_equal(stringsim(x,y,method="osa", useBytes=TRUE ), 1-2/3)
  expect_equal(stringsim(x,y,method="lv", useBytes=FALSE),  1-1/2)
  expect_equal(stringsim(x,y,method="lv", useBytes=TRUE ),  1-2/3)
  expect_equal(stringsim(x,y,method="dl", useBytes=FALSE),  1-1/2)
  expect_equal(stringsim(x,y,method="dl", useBytes=TRUE ),  1-2/3)
  expect_equal(stringsim(x,y,method="hamming", useBytes=FALSE), 1-1/2)
  expect_equal(stringsim(x,y,method="hamming", useBytes=TRUE ), 1-1)
  expect_equal(stringsim(x,y,method="lcs", useBytes=FALSE),     1-1/2)
  expect_equal(stringsim(x,y,method="lcs", useBytes=TRUE ),     1-3/5)
  expect_equal(stringsim(x,y,method="qgram", q=1, useBytes=FALSE), 1-1/2)
  expect_equal(stringsim(x,y,method="qgram", q=1, useBytes=TRUE ), 1-3/5)
  expect_equal(stringsim(x,y,method="cosine", q=1, useBytes=FALSE), 1-1/2)
  expect_equal(stringsim(x,y,method="cosine", q=1, useBytes=TRUE ), 1-(1-1/sqrt(6)))
  expect_equal(stringsim(x,y,method="jaccard", q=1, useBytes=FALSE), 1-2/3)
  expect_equal(stringsim(x,y,method="jaccard", q=1, useBytes=TRUE ), 1-3/4)
  expect_equal(stringsim(x,y,method="jw", useBytes=FALSE), 1-1/3)
  expect_equal(stringsim(x,y,method="jw", useBytes=TRUE ), (1/2 + 1/3 +1)/3)

## seq_sim

## elementary seq_sim test
  expect_equal(
    seq_sim(list(1:3,2:4),list(1:3))
     , stringsim(c("abc","bcd"),"abc") 
  )






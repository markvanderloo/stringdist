
library(testthat)

### -------------------------------------------------------------
context("Phonetic")
test_that("Soundex",{

  testset <- "name;code
Robert;R163
rupert;R163
Rubin;R150
Ashcraft;A261
asHCroft;A261
Tymczak;T522
Pfister;P236
gutierrez;G362
Jackson;J250
washington;W252
Lee;L000
NA;NA"
  testset <- read.csv2(textConnection(testset), stringsAsFactors=FALSE)
  expect_that(phonetic(testset$name,"soundex"), equals(testset$code))

}) 




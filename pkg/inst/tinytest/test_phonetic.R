options(sd_num_thread=2)

### -------------------------------------------------------------

## Soundex

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
  expect_equal(phonetic(testset$name,"soundex"), testset$code)
  expect_equal(phonetic(testset$name,"soundex",useBytes=TRUE), testset$code)
  expect_warning(phonetic(paste0('Mot',intToUtf8(0x00F6),'rhead')))  


## soundex handles encoding
  ouml <- intToUtf8("0x00F6")
  # non-ascii within string
  expect_warning(phonetic(paste0("Mot",ouml,"rhead"),method='soundex'))
  # non-ascii at beginning of string
  expect_warning(phonetic(paste0(ouml,"zzy"),method='soundex'))
  # non-printable in string (carriage return)
  cr <- "\r"
  expect_warning(phonetic(paste0(cr,"hello"),method='soundex'))





options(sd_num_thread=1L)

# tests against cases that used to segfault when we did not check
# NULL cases.
expect_error(afind("a","b",nthread=1:4))
expect_error(afind("a","b",nthread="foo"))
expect_error(afind("a","b",nthread=integer(0)))
expect_error(afind("a","b",nthread=NULL))





texts = c("When I grow up, I want to be"
  , "one of the harversters of the sea"
  , "I think before my days are gone"
  , "I want to be a fisherman")

patterns = c("fish", "gone","to be")

out <- afind(texts, patterns, method="osa")

location <- matrix(c(
   1,    1,  24,
   6,    1,  28,
   1,   28,   6,
  16,    3,   8),
  nrow=4, byrow=TRUE)


distance <- matrix(c(
   4,    3,    0,
   2,    2,    3,
   3,    0,    2,
   0,    3,    0),
   nrow=4, byrow=TRUE)

match <- matrix(c(
  "When", "When", "to be",
  "f th", "one ", "he se",
  "I th", "gone", "nk be",
  "fish", "want", "to be"),
   nrow=4, byrow=TRUE)


expect_equal(out$location, location)
expect_equal(out$distance, distance)
expect_equal(out$match,    match)

# test paralellization

out1 <- afind(texts, patterns, method="osa", nthread=2L)
expect_identical(out, out1)

# test 'value' option
out2 <- afind(texts, patterns, value=FALSE)
expect_equal(length(out2), 2)


# test grep/grepl equivalents 'grab', 'grabl'

expect_equal(grab(texts, "harvester", maxDist=2), 2)
expect_equal(grab(texts, "harvester", value=TRUE, maxDist=2), "harverste")
expect_equal(grabl(texts, "harvester", maxDist=2)
            , c(FALSE,TRUE,FALSE,FALSE))

expect_equal(extract(texts, "harvester", maxDist=2)
            , matrix(c(NA, "harverste",NA,NA),nrow=4) )

## Test running_cosine
pattern <- c("phish", "want to")

expect_identical(
    afind(texts, pattern, method="cosine",         q=3)
  , afind(texts, pattern, method="running_cosine", q=3)
)


## test whether the correct positions are returned for all methods.

methods = names(stringdist:::METHODS)
methods = methods[!methods %in% c("soundex","hamming")]
text <- "If you squeeze my lizzard, I put my snake on you."
pattern <- "lizard"

for ( method in methods ){
  expect_equal(afind(text, pattern, method=method, q=3, p=0.1)$location[1,1], 19, info=method)
}

## test the usual edge cases

# notice: window size = 0.
expect_equal(afind("foo","")$distance[1], 0)

expect_equal(afind("foo",NA)$distance[1], NA_real_)
expect_equal(afind("foo",NA)$location[1], NA_integer_)
expect_equal(afind("foo",NA)$match[1], NA_character_)

expect_equal(afind(NA,"foo")$distance[1], NA_real_)
expect_equal(afind(NA,"foo")$location[1], NA_integer_)
expect_equal(afind(NA,"foo")$match[1], NA_character_)

expect_equal(afind("","foo")$distance[1], 3)
expect_equal(afind("","foo")$location[1], 1)
expect_equal(afind("","foo")$match[1], "")

expect_equal(grab("foo", ""), 1L)
expect_equal(grabl("foo",""), TRUE)
expect_equal(grab("foo",NA), integer(0))

# note that 'grepl' gives FALSE in this case (which is inconsistent with
# grepl(NA, NA), grepl(NA, "foo").
expect_equal(grabl("foo",NA), NA)















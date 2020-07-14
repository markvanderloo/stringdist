options(sd_num_thread=1L)

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

# test option

out2 <- afind(texts, patterns, value=FALSE)
expect_equal(length(out2), 2)


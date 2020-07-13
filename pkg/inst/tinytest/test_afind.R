


x <- c(
  "when I grow up I want to be, one of the harverster of the sea"
, "I think before my days are gone, I want to be a fisherman"
)

out <- afind(x, "gone")
expect_equal(out$location, matrix(c(29,28),nrow=2))
expect_equal(out$match, matrix(c(" one", "gone"),nrow=2))




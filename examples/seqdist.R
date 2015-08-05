# Distances between lists of integer vectors. Note the postfix 'L' to force
# integer storage.
a <- list(c(102L, 107L))     # fu
b <- list(c(102L,111L,111L)) # foo
seqdist(a,b)

# also works for stringdistmatrix
a <- lapply(c("foo","bar","baz"),utf8ToInt)
seqdistmatrix(a)


# Distances between lists of integer vectors. Note the postfix 'L' to force 
# integer storage. The shorter argument is recycled over (\code{a})
a <- list(c(102L, 107L))                        # fu
b <- list(c(102L,111L,111L),c(102L,111L,111L))  # foo, fo
seq_dist(a,b)

# translate strings to a list of integer sequences 
a <- lapply(c("foo","bar","baz"),utf8ToInt)
seq_distmatrix(a)

# Note how missing values are treated. NA's as part of the sequence are treated 
# as an integer (the representation of \code{NA_integer_}).
a <- list(NA_integer_,c(102L, 107L))
b <- list(c(102L,111L,111L),c(102L,111L,NA_integer_))  
seq_dist(a,b)

\dontrun{
# Distance between sentences based on word order. Note: words must match exactly or they
# are treated as completely different.
#
# For this example you need to have the 'hashr' package installed.
a <- "Mary had a little lamb"
a.words <- strsplit(a,"[[:blank:]]+")
a.int <- hashr::hash(a.words)
b <- c("a little lamb had Mary",
           "had Mary a little lamb")
b.int <- hashr::hash(strsplit(b,"[[:blank:]]+"))
seq_dist(a.int,b.int)
}


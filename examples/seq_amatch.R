
x <- list(1:3,c(3:1),c(1L,3L,4L))
table <- list(
  c(5L,3L,1L,2L)
  ,1:4
)
seq_amatch(x,table,maxDist=2)

# behaviour with missings
seq_amatch(list(c(1L,NA_integer_,3L),NA_integer_), list(1:3),maxDist=1)


\dontrun{
# Match sentences based on word order. Note: words must match exactly or they
# are treated as completely different.
#
# For this example you need to have the 'hashr' package installed.
x <- "Mary had a little lamb"
x.words <- strsplit(x,"[[:blank:]]+")
x.int <- hashr::hash(x.words)
table <- c("a little lamb had Mary",
           "had Mary a little lamb")
table.int <- hashr::hash(strsplit(table,"[[:blank:]]+"))
seq_amatch(x.int,table.int,maxDist=3)
}

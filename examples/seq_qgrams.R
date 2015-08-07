
# compare the 2-gram overlap between sequences 1:3 and 2:4
seq_qgrams(x = 1:3, y=2:4,q=2)

# behavior when NA's are present.
seq_qgrams(1:3,c(1,NA,2),NA_integer_)

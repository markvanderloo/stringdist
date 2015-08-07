L1 <- list(1:3,2:4)
L2 <- list(1:3)
seq_sim(L1,L2,method="osa")

# note how missing values are handled (L2 is recycled over L1)
L1 <- list(c(1L,NA_integer_,3L),2:4,NA_integer_)
L2 <- list(1:3)
seq_sim(L1,L2)


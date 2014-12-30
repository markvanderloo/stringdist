# smoke tests.
set.seed(1864)

dyn.load("../pkg/src/stringdist.so")
for (f in dir("../pkg/R/",full.names=TRUE)) d <- source(f)
options(sd_num_thread=4)

# printable ascii
pascii <- sapply(33:126,intToUtf8)
# the ascii set
ascii <- sapply(0:127,intToUtf8)
# the utf8 set
utf8 <- sapply(0:0x10FFFF,intToUtf8)
utf8 <- utf8[stringi::stri_enc_isutf8(utf8)]


# methods
edit <- c("osa", "lv", "dl", "lcs", "hamming") 
qgrm <- c("qgram","cosine", "jaccard")
heur <- c("jw")
phon <- c("soundex")

useBytes <- c(TRUE,FALSE)
nthread=1:4
q <- 0:3
p <- c(0,0.01,0.25)

tests <- list(
    test_edit = expand.grid(method = edit, useBytes = useBytes, nthread=1:4,stringsAsFactors = FALSE)
  , test_qgrm = expand.grid(method = qgrm, useBytes = useBytes, nthread=1:4,q = q, stringsAsFactors = FALSE)
  , test_heur = expand.grid(method = heur, useBytes = useBytes, nthread=1:4,p = p, stringsAsFactors = FALSE)
  , test_phon = expand.grid(method = phon, useBytes = useBytes, nthread=1:4,stringsAsFactors=FALSE)
)

rand_string <- function(n, len, alphabet=pascii){
  stopifnot( length(len) == 1 || length(len) == n )
  if (length(len)==1){
    sapply(1:n,function(i) paste0(sample(alphabet, len, replace=TRUE),collapse=""))
  } else {
    sapply(1:n,function(i) paste0(sample(alphabet, len[i], replace=TRUE),collapse=""))    
  }  
}

N1 <- 250
N2 <- N1
len1 <- sample(0:50,N1,replace=TRUE)
len2 <- sample(0:50,N2,replace=TRUE)
# smoke test for edit-based distances
str1 <- rand_string(N1,len1,pascii)
str2 <- c(rand_string(N2,len2,pascii)
          , rand_string(N2,len2,ascii)
          , rand_string(N2,len2,utf8))

for (tst in tests ){
  for ( i in 1:nrow(tst) ){
    args <- c(list(a=str1,b=str2),as.list(tst[i,]))
    tryCatch(do.call(stringdist,args), error=function(e){
      cat(sprintf("Failure with message\n%s\n",e$message))
      print(tst[i,])
    }, warning=function(w){
      cat(sprintf("Warning with message\n%s\n",w$message))
      print(tst[i,])      
    })
  }
}

for (tst in tests ){
  for ( i in 1:nrow(tst) ){
    args <- c(list(x=str1,table=str2),as.list(tst[i,]))
    tryCatch(do.call(amatch,args), error=function(e){
      cat(sprintf("Failure with message\n%s\n",e$message))
      print(tst[i,])
    }, warning=function(w){
      cat(sprintf("Warning with message\n%s\n",w$message))
      print(tst[i,])      
    })
  }
}



#tst <- tests[[4]]
#args <- c(list(a=str1[1],b=str2[1]),as.list(tst[1,]))
#do.call(stringdist,args)









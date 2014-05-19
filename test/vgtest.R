
dyn.load("../pkg/src/stringdist.so")
for ( f in dir("../pkg/R",full.names=TRUE) ) source(f)

 x <- paste0('Mot',intToUtf8(0x00F6),'rhead')
 y <- c('bastard','Motorhead')
 jwdist <- round(1-(1/3)*(8/9 + 8/10 + 1),3)

#amatch(x, y, method='dl', maxDist=2, useBytes=TRUE)
#stringdist(x[1], y[2], method='dl', useBytes=FALSE)
#stringdist(x[1], y[2], method='dl', useBytes=TRUE)
#stringdist("b","a" , method='dl',weight=c(1,1,0.5,1))
stringdist("b","a" , method='dl',weight=c(1,1,1,1))
#stringdist(x, y, method='osa', useBytes=TRUE)



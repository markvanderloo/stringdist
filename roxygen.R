library(roxygen2)
library(pkgload)
options(error=traceback)
unlink( 'pkg/man', TRUE)

#setwd('pkg')
load_all('./pkg')
document('./pkg')
#roxygenize( '.'
#          , roxygen.dir='.'
#          , copy.package=FALSE
#          , unlink.target=TRUE
#		    )


if (length(list.files('inst/doc')) == 0){
   unlink( 'inst/doc', TRUE)   
}

#' A package for string distance calculation
#' 
#' @name stringdist-package
#' @docType package
#' @useDynLib stringdist
{}

#' Compute distance between strings
#'
#' @section Details:
#' Currently, the following distance measures are supported:
#' \tabular{ll}{
#'    \code{osa} \tab Optimal string aligment, (restricted Damerau-Levenshtein distance).\cr
#'    \code{lv} \tab Levenshtein distance.\cr
#'    \code{dl} \tab Full Damerau-Levenshtein distance.\cr
#'    \code{h}  \tab Hamming distance (\code{a} and \code{b} must be of equal size).
#' }
#' The Hamming distance counts the number of character substitutions that turns \code{a} into \code{b}. The Levenshtein distance
#' allows deletions, insertions and substitutions. The Optimal string alignment distance also allows transpositions, but each substring
#' may be edited only once, so a character cannot be transposed twice. The Damerau-Levensthein distance alows multiple transpositions.
#'
#' @section Encoding issues:
#' Input strings are re-encoded to \code{utf8} an then to \code{integer}
#' vectors prior to the distance calculation (which works on unsigned ints). 
#' Unfortunately this double conversion is necessary as it seems the only way to
#' reliably convert (possibly multibyte) characters to integers on all systems
#' supported by \code{R}.
#' (\code{R}'s native \code{\link[utils]{adist}} function does this as well). See \code{\link[base]{Encoding}} for further details.
#'
#'
#' @param a \code{character} vector
#' @param b \code{character} vector
#' @param method Method for distance calculation (see details)
#' @param weight The penalty for deletion, insertion, substitution and transposition (where applicable).
#' @param maxDist Maximum string distance before calculation is stopped, \code{maxDist=0} means calculation goes on untill the distance is computed.
#' @return For \code{stringdist},  a vector with string distances of size \code{max(length(a),length(b))}.
#'  For \code{stringdistmatrix}, a \code{length(a)xlength(b)} \code{matrix}. The returned distance is \code{-1} when \code{maxDist} is exceeded
#'  and \code{NA} if any of \code{a} or \code{b} is \code{NA}.
#' @export
stringdist <- function(a, b, method=c("osa","lv","dl","h"), weight=c(d=1,i=1,s=1,t=1), maxDist=0){
   a <- as.character(a)
   b <- as.character(b)
   if (length(a) == 0 || length(b) == 0){ 
      return(numeric(0))
   }
   method <- match.arg(method)
   a <- lapply(enc2utf8(a),utf8ToInt)
   b <- lapply(enc2utf8(b),utf8ToInt)
   do_dist(a,b,method,weight,maxDist)
}


do_dist <- function(a,b,method,weight,maxDist){
   switch(method,
      osa = .Call('R_osa', a, b, as.double(weight), as.double(maxDist)),
      lv  = .Call('R_lv' , a, b, as.double(weight), as.double(maxDist)),
      dl  = .Call('R_dl' , a, b, as.double(weight), as.double(maxDist)),
      h   = .Call('R_hm' , a, b, as.integer(maxDist))
   )
}


#' @param ncores number of cores to use. Parallelisation is over \code{b}, so the speed gain by parallelisation is highest when \code{b} is shorter than \code{a}.
#' @rdname stringdist
#' @export
stringdistmatrix <- function(a, b, method=c("osa","lv","dl","h"), weight=c(d=1,i=1,s=1,t=1), maxDist=0,ncores=1){
   a <- as.character(a)
   b <- as.character(b)
   if (length(a) == 0 || length(b) == 0){ 
      return(numeric(0))
   }
   method <- match.arg(method)
   a <- lapply(enc2utf8(a),utf8ToInt)
   b <- lapply(enc2utf8(b),utf8ToInt)
   if (ncores==1){
     x <- sapply(b,do_dist,a,method,weight,maxDist)
   } else {
      cl <- makeCluster(ncores)
         x <- parSapply(cl, b,do_dist,a,method,weight,maxDist)
      stopCluster()
   }
   x
}




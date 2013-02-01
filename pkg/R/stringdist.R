#' A package for string distance calculation
#' 
#' @name stringdist-package
#' @docType package
#' @useDynLib stringdist
{}

#' Compute distance between two strings
#'
#' @section Details:
#' Computes pairwise distances between two \code{character} vectors. The shortest of \code{a} and \code{b} is recycled.
#'
#' @param a \code{character} vector
#' @param b \code{character} vector
#' @param method Method for distance calculation (see details)
#' @param weight The penalty for deletion, insertion, substitution and transposition (where applicable).
#' @param maxDist Maximum string distance before calculation is stopped, \code{maxDist=0} means calculation goes on untill the distance is computed.
#' @return A vector with string distances of size \code{max(length(a),length(b))}
#' @export
stringdist <- function(a, b, method=c("osa","dl","h"), weight=c(d=1,i=1,s=1,t=1), maxDist=0){
   a <- as.character(a)
   b <- as.character(b)
   if (length(a) == 0 || length(b) == 0){ 
      return(numeric(0))
   }
   method <- match.arg(method)
   a <- lapply(enc2utf8(a),utf8ToInt)
   b <- lapply(enc2utf8(b),utf8ToInt)
   switch(method,
      osa = .Call('R_osa', a, b, as.double(weight), as.double(maxDist)),
      dl  = .Call('R_dl' , a, b, as.double(weight), as.double(maxDist)),
      h   = .Call('R_hm' , a, b, as.integer(maxDist))
   )
}



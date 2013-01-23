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
stringdist <- function(a, b, method=c("osa","dl"), weight=c(d=1,i=1,s=1,t=1), maxDist=0){
   stopifnot(is.character(a), is.character(b))
   if (length(a) == 0 || length(b) == 0){ 
      return(numeric(0))
   }
   method <- match.arg(method)
   na <- nchar(a)
   nb <- nchar(b)
   switch(method,
      osa = .Call('R_osa', a, b, na, nb, as.double(weight), as.double(maxDist)),
      dl  = .Call('R_dl' , a, b, na, nb, as.double(weight), as.double(maxDist), max(max(na),max(nb)))
   )
}



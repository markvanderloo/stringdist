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
#' @return A vector with string distances of size \code{max(length(a),length(b))}
#' @export
stringdist <- function(a, b, method="osa", weight=c(d=1,i=1,s=1,t=1)){
   stopifnot(is.character(a), is.character(b))

   .Call('R_osa',
      a,
      b,
      nchar(a),
      nchar(b),
      as.double(weight)
   )
}



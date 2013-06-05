#' Get a table of qgram counts from a character vector.
#'
#' @section Details:
#' The input is converted to \code{character}. Each element is converted to \code{integer} via \code{utf8} 
#' as in \code{\link{stringdist}}, prior to passing the data to the underlying routine. The \code{names} of the output 
#' table (i.e. the qgrams) are encoded in \code{utf8}.
#'
#'
#'
#' @param x character vector
#' @param q size of q-gram, must be non-negative.
#'
#' @return An object of class \code{table}
#'
#' @seealso \code{\link{stringdist}}. 
#'
#' @example ../examples/qgrams.R
#' @export
qgrams <- function(x, q ){
  x <- as.character(x)
  q <- as.integer(q)

  v <- .Call("R_get_qgrams", char2int(x), q)
  if ( is.null(v) || length(v) == 0 ) return( table(integer(0)) )

  A <- array(attr(v,"qgrams"),dim=c(q,length(v)))
  qgrams <- apply(A,2,intToUtf8)
  attr(v,"qgrams") <- NULL
  names(v) <- qgrams
  as.table(v)
}



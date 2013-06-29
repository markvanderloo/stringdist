#' Get a table of qgram counts from a character vector.
#'
#' @section Details:
#' The input is converted to \code{character}. Each element is converted to \code{utf8} and than to \code{integer}  
#' as in \code{\link{stringdist}}. Next,the data is passed to the underlying routine. row names of the output 
#' array (i.e. the qgrams) are encoded in \code{utf8}.
#' 
#' Strings with less than \code{q} characters and elements containing \code{NA} are skipped. Using \code{q=0} 
#' therefore counts the number of empty strings \code{""} occuring in each argument.
#' 
#' @param ... any number of (named) arguments, that will be coerced to character with \code{as.character}.
#' @param q size of q-gram, must be non-negative.
#'
#' @return An table with $q$-gram counts. Detected q-grams are column names and the argument names as row names.
#' If no argument names were provided, they will be generated.
#'
#' @seealso \code{\link{stringdist}}. 
#'
#' @example ../examples/qgrams.R
#' @export
qgrams <- function(...,q=1L){
  q <- as.integer(q)
  L <- lapply(list(...), as.character)

  if ( is.null(names(L)) ){
    names(L) <- paste("V",1:length(L),sep="")
  } else {
    I <- names(L) == "";
    names(L)[I] = paste("V",which(I),sep="")
  }
  L <- lapply(L,char2int)
  
  v <- .Call("R_get_qgrams",L,as.integer(q))
  # get qgram labels
  A <- array(attr(v,"qgrams"),dim=c(1,length(v)))
  if ( q == 0 ){
    qgrams = ""
  } else {
    Q <- attr(v,"qgrams")
    attr(v,"qgrams") <- NULL
    A <- array(Q,dim=c(q, length(v)/length(L)))
    qgrams = apply(A,2,intToUtf8)
  }
  array(v,
    dim=c(length(L), length(qgrams)),
    dimnames = list(
    names(L),
    qgrams
   ) 
  )
}



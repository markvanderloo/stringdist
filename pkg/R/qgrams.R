#' Get a table of qgram counts from one or more character vectors.
#'
#' @section Details:
#' The input is converted to \code{character}. If \code{useBytes=TRUE}, each element is 
#' converted to \code{utf8} and then to \code{integer} as in \code{\link{stringdist}}. 
#' Next,the data is passed to the underlying routine. row names of the output 
#' array (i.e. the qgrams) are encoded in \code{utf8}.
#' 
#' Strings with less than \code{q} characters and elements containing \code{NA} are skipped. Using \code{q=0} 
#' therefore counts the number of empty strings \code{""} occuring in each argument.
#' 
#' @param ... any number of (named) arguments, that will be coerced to character with \code{as.character}.
#' @param q size of q-gram, must be non-negative.
#' @param useBytes Determine byte-wise qgrams. \code{useBytes=TRUE} is faster but may yield different
#' 	results depending on character encoding. For \code{ASCII} it is identical. See also \code{\link{stringdist}} under Encoding issues.
#' @param useNames Add q-grams as column names. If \code{useBytes=useNames=TRUE}, the q-byte sequences are represented as 2 hexadecimal numbers
#'   per byte, separated by a vertical bar (\code{|}).
#' @param .list Will be concatenated with the \code{...} argument(s). Usefull for adding character vectors named \code{'q'} or \code{'useNames'}.
#' @return A table with \eqn{q}-gram counts. Detected \eqn{q}-grams are column names and the argument names as row names.
#' If no argument names were provided, they will be generated.
#'
#' @seealso \code{\link{stringdist}}, \code{\link{amatch}}
#'
#' @example ../examples/qgrams.R
#' @export
qgrams <- function(..., .list=NULL,q=1L,useBytes=FALSE, useNames=!useBytes){
  q <- as.integer(q)

  if (!is.null(.list) && length(.list) == 0) .list=NULL
  L <- lapply(c(list(...),.list), as.character)
  if (length(L) == 0) return(array(dim=c(0,0)))

  if ( is.null(names(L)) ){
    names(L) <- paste("V",1:length(L),sep="")
  } else {
    I <- names(L) == "";
    names(L)[I] = paste("V",which(I),sep="")
  }
  L <- lapply(L,char2int)
  
  v <- .Call("R_get_qgrams",L,as.integer(q))
  
  nqgrams <- length(v)/length(L)
  qgrams <- NULL
  if (useNames){  
    if ( q == 0 ){
      qgrams = ""
    } else {
      Q <- attr(v,"qgrams")
      attr(v,"qgrams") <- NULL
      A <- array(Q,dim=c(q, nqgrams))
      qgrams = if( useBytes ){
            apply(A,2,function(x) paste(as.raw(x),collapse="|")) 
          } else {
            apply(A,2,intToUtf8)
          }
    }
  }
  array(v,
    dim=c(length(L), nqgrams),
    dimnames = list(
    names(L),
    qgrams
   ) 
  )
}


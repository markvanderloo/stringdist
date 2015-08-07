#' Get a table of qgram counts from one or more character vectors.
#'
#' @section Details:
#' The input is converted to \code{character}. If \code{useBytes=TRUE}, each element is 
#' converted to \code{utf8} and then to \code{integer} as in \code{\link{stringdist}}. 
#' Next,the data is passed to the underlying routine.
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
#' @param .list Will be concatenated with the \code{...} argument(s). Useful for adding character vectors named \code{'q'} or \code{'useNames'}.
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
  L <- setnames(L)
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
            enc2native(apply(A,2,intToUtf8))
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


setnames <- function(x){
  list_names <- names(x)
  generic_names <- paste0("V",seq_along(x))
  if (is.null(list_names)){
    return(setNames(x,generic_names))
  }
  I <- list_names == ""
  names(x)[I] <- generic_names[I]
  x
}

#' Get a table of qgram counts for integer sequences
#' 
#' 
#' @param ... Any number of (named) arguments that will be coerced with \code{as.integer}
#' @param .list Will be concatenated with the \code{...} argument(s). Useful for adding integer vectors named 'q'.
#' @param q The size of q-gramming.
#' 
#' @return 
#' A \code{matrix} containing q-gram profiles. Columns 1 to \code{q} contain the
#' encountered q-grams. The ensuing (named) columns contain the q-gram counts
#' per vector. Run the example for a simple overview.
#' 
#' Missing values in integer sequences are treated as any other number.  
#' 
#' @example ../examples/seq_qgrams.R
#' 
#' @seealso \code{\link{seq_dist}}, \code{\link{seq_amatch}}
#' 
#' @export 
#' 
seq_qgrams <- function(...,.list=NULL,q=1L){
  L <- lapply(c(list(...),.list),function(x) list(as.integer(x)))
  if (length(L) == 0) return(array(dim=c(0,0)))
  L <- setnames(L)
  v <- .Call("R_get_qgrams",L,as.integer(q))
  Q <- attr(v,"qgrams")
  nqgrams <- length(v)/length(L)
  Q <- t(array(Q,dim=c(q,nqgrams),dimnames=list(paste0("q",1:q),NULL)))
  A <- matrix(v, nrow=nqgrams, ncol=length(L),byrow=TRUE,dimnames=list(NULL,paste0("n.",names(L))))
  cbind(Q,A)
}

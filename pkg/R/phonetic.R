
#' Phonetic algorithms
#' 
#' Translate strings to phonetic codes. Similar sounding strings should get 
#' similar or equal codes. 
#'
#' @param x a character vector whose elements are phonetically encoded. 
#' @param method name of the algorithm used. The default is \code{"soundex"}.
#' @param useBytes Perform byte-wise comparison. \code{useBytes=TRUE} is faster 
#'    but may yield different results depending on character encoding. For more
#'    information see the documentation of \code{\link{stringdist}}.
#'
#' @details
#' Currently, only the soundex algorithm is implemented. Note that soundex coding
#' is only meaningful for characters in the ranges a-z and A-Z. Specifically, soundex
#' coding of strings containing non-ascii characters may be system-dependent and should 
#' not be trusted. 
#' 
#' 
#' @return
#' The returns value depends on the method used. However, all currently 
#' implemented methods return a character vector of the same length of the input
#' vector. 
#'
#' @section Citation:
#' If you would like to cite this package, please cite the R-journal paper: 
#' \itemize{
#' \item{M.P.J. van der Loo (2014). The \code{stringdist} package for approximate string matching. 
#'  R Journal 6 (accepted for publication)}
#' }
#' Or use \code{citation('stringdist')} to get a bibtex item.
#' 
#' @references
#' \itemize{
#' \item{The Soudex algorithm implemented is the algorithm used by the 
#'   \href{http://www.archives.gov/research/census/soundex.html}{National Archives}. 
#'   This algorithm differs slightly from the original algorithm patented by R.C. Russell 
#'   (US patents 1261167 (1918) and 1435663 (1922)). 
#' }
#' }
#'
#' @example ../examples/phonetic.R
#'
#' @export
phonetic <- function(x, method = c("soundex"), useBytes = FALSE) {
  x <- as.character(x)
  method <- match.arg(method)
  stopifnot(is.logical(useBytes))
  if (!useBytes) x <- char2int(x)
  if (method == "soundex") {
    r <- .Call("R_soundex", x)
    if (!useBytes) int2char(r) else r
  } 
}

int2char <- function(x) {
  enc2native(sapply(x, intToUtf8))
}


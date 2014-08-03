
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
#' Currently, only the soundex algorithm is implemented. 
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
#'
#' @export
phonetic <- function(x, method = c("soundex"), useBytes = FALSE) {
  x <- as.character(x)
  method <- match.arg(method)
  stopifnot(is.logical(useBytes))
  if (method == "soundex") {
    .Call("R_soundex", x)
  }
}


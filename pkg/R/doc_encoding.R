#' @title 
#' String metrics in \pkg{stringdist}
#' 
#' @description
#' This page gives an overview of encoding handling in \pkg{stringst}. 
#'
#' 
#' @section Encoding in \pkg{stringdist}:
#' 
#' All character strings are stored as a sequence of bytes. An encoding
#' system relates a byte, or a short sequence of bytes to a symbol. Over the years, many 
#' encoding systems have been developed, and not all OS's and softwares use the same encoding 
#' as default. Similarly, depending on the system R is running on, R may use a
#' different encoding for storing strings internally.
#'
#' The \pkg{stringdist} package is designed so users in principle need not
#' worry about this. Strings are converted to \code{UTF-32} (unsigned integer)
#' by default prior to any further computation. This means that results are
#' encoding-independent and that strings are interpreted as a sequence of
#' symbols, not as a sequence of pure bytes. In functions where this is
#' relevant, this may be switched by setting the \code{useBytes} option to
#' \code{TRUE}. However, keep in mind that results will then likely depend on the
#' system R is running on, except when your strings are pure ASCII.
#' Also, for multi-byte encodings, results for byte-wise computations
#' will usually differ from results using encoded computations.
#'
#' Prior to \pkg{stringdist} version 0.9, setting \code{useBytes=TRUE} could 
#' give a significant performance enhancement. Since version 0.9, translation
#' to integer is done by C code internal to \pkg{stringdist} and the difference in
#' performance is now negligible.
#'
#' @section Unicode normalisation:
#' In \code{utf-8}, the same (accented) character may be represented as several byte sequences. For example, an u-umlaut
#' can be represented with a single byte code or as a byte code representing \code{'u'} followed by a modifier byte code
#' that adds the umlaut. The \href{http://cran.r-project.org/package=stringi}{stringi} package 
#' of Gagolevski and Tartanus offers unicode normalisation tools. 
#'
#' @section Some tips on character encoding and transliteration:
#' Some algorithms (like soundex) are defined only on the printable ASCII character set. This excludes any character
#' with accents for example. Translating accented characters to the non-accented ones is a form of transliteration. On
#' many systems running R (but not all!) you can achieve this with 
#' 
#' \code{iconv(x,to="ASCII//TRANSLIT")}, 
#' 
#' where \code{x} is your character vector. See the documentation of \code{\link[base]{iconv}} for details.
#'
#' The \code{stringi} package (Gagolewski and Tartanus) should work on any system. The command 
#' \code{stringi::stri_trans_general(x,"Latin-ASCII")} transliterates character vector \code{x} to ASCII.
#'
#' @references
#' \itemize{
#'  \item{The help page of \code{\link[base]{Encoding}}} describes how R handles encoding.
#'  \item{The help page of \code{\link[base]{iconv}} has a good overview of base R's 
#'       encoding conversion options. The capabilities of \code{iconv} depend on the system R is running on.
#'       The \pkg{stringi} package offers platform-independent encoding and normalization tools.}
#' }
#'
#' @seealso
#' \itemize{
#' \item{Functions using re-encoding: \code{\link{stringdist}}, \code{\link{stringdistmatrix}}, \code{\link{amatch}}, \code{\link{ain}}, \code{\link{qgrams}}}
#' \item{Encoding related: \code{\link{printable_ascii}}}
#' }
#' @name stringdist-encoding
#' @rdname stringdist-encoding
{}


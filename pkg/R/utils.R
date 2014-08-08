
RECYCLEWARNING <- NULL


.onLoad <- function(libname, pkgname){
  RECYCLEWARNING <<- gettext(tryCatch( (1:2)+(1:3),warning=function(w) w$message ))
}

#' Detect the presence of non-printable or non-ascii characters
#' 
#' @param x a \code{character} vector
#'
#' @details
#' Printable ASCII characters consist of space, \code{A-Z}, \code{a-z}, \code{0-9} and the characters
#'
#' \code{! "" # $ \% & ' ( ) * + , . / : ; < = > ? @@ [ ] \\ ^ _ ` { | } ~ -} 
#'
#' Note that this excludes tab (as it is a control character).
#'
#' @section Some tips on character encoding and transliteration:
#' Some algorithms (like soundex) are defined only on the printable ASCII character set. This excludes any character
#' with accents for example. Translating accented characters to the non-accented ones is a form of transliteration. On
#' most systems running R you can achieve this with 
#' 
#' \code{iconv(x,"ASCII//TRANSLIT")}, 
#' 
#' where \code{x} is your character vector. See the documentation of \code{\link[base]{iconv}} for details.
#'
#' The \code{stringi} package (Gagolewski and Tartanus) should work on any system. The command 
#' \code{stringi::stri_trans_general(x,"Latin-ASCII")} transliterates character vector \code{x} to ASCII.
#' 
#' @example ../examples/printable_ascii.R
#'
#'
#' @return A \code{logical} indicating which elements consist solely of printable ASCII characters.
#' @export 
printable_ascii <- function(x){ 
  # To portably detect ASCII characters, we need to specify them literally. Hence this monster of a character class. See ?regexp.

  # notes:
  # - caret (^) at the beginning negates what comes after, the caret in the middle is the actual character.
  # - the closing square bracket ] needs to be specified first
  # - double quote " and backslash are escaped
  # - the dash "-" is specified at the end since it would indicate a range otherwise
  # - see ? regexp.
  charclass <- paste0("[^]"
    , "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    , "0123456789"
    , " !\"#$%&'()*+,./:;<=>?@[\\^_`{|}~-"
    , "]"
  )  
  !grepl(charclass,x)
}




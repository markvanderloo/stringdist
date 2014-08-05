
RECYCLEWARNING <- NULL


.onLoad <- function(libname, pkgname){
  RECYCLEWARNING <<- gettext(tryCatch( (1:2)+(1:3),warning=function(w) w$message ))
}

#' Check if a character consists of printable ASCII characters.
#' 
#' @param x a \code{character} vector
#'
#' @details
#' Printable ASCII characters consist of \code{A-Z}, \code{a-z}, \code{0-9} and the punctuation characters
#'
#' \code{^ ( ) ! \" \\  # $ \%  & ' + * , . / : ; < = > ? @@} 
#'
#' @return A \code{logical} indicating which elements consist solely of printable ASCII characters.
#' @export 
printable_ascii <- function(x){ 
  # To portably detect ASCII characters, we need to specify them literally. Hence this monster of a character class. See ?regexp.
  # Also note: ^ in the beginning negates,^ further on is a character. Literal [ must be specified at the front of the list
  # and literal dash '-' must be specified last. 
  !grepl("[^][ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789^()!\"#$%&'+*,./:;<=>?@\\_`{}|~-]",x)
}



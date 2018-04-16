.onUnload <- function (libpath) {
  library.dynam.unload("stringdist", libpath)
}

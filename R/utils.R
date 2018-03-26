fprintf <- function(fmt, ..., file = "", append = FALSE) {
  mystr <- sprintf(fmt, ...)
  cat(mystr, file = file, append = append)
  invisible(nchar(mystr))
}

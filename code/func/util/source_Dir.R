## If you want to source() a bunch of files, something like
## the following may be useful:
sourceDir <- function(path, trace = TRUE, ...) {
  nm_lst <- 
    list.files(path, pattern = "[.][R]$", recursive = TRUE)
  for (nm in nm_lst) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
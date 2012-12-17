### =========================================================================
### Utilities
### -------------------------------------------------------------------------

## May push up to gmapR if it proves more generally useful

flankingCycleBreaks <- function(read_length, width = 10L) {
  if (is.na(read_length))
    return(NULL)
  if (read_length < 1)
    stop("'read_length' must be >= 1 or NA")
  if (width < 0)
    stop("'width' must be non-negative")
  
  as.integer(c(0L, width, read_length - width, read_length))
}

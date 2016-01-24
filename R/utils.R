### =========================================================================
### Utilities
### -------------------------------------------------------------------------

## May push up to gmapR if it proves more generally useful

flankingCycleBreaks <- function(read_length, width = 10L) {
  read_length <- as.integer(read_length)
  if (is.na(read_length))
    return(NULL)
  if (read_length < 1)
    stop("'read_length' must be >= 1 or NA")
  if (width < 0)
    stop("'width' must be non-negative")
  
  as.integer(c(0L, width, read_length - width, read_length))
}

installed <- function(x) {
  !identical(suppressWarnings(packageDescription(x)), NA)
}

chunkRange <- function(which, n) {
  chunks <- breakInChunks(width(which), ceiling(width(which) / n))
  which <- GRanges(seqnames(which),
                   IRanges(start(which) + start(chunks) - 1L,
                           width = width(chunks)))

}

## some internal compatibility wrappers for older versions
rawTotalDepth <- function(x) {
    if (!is.null(x$raw.count.total))
        x$raw.count.total
    else x$count.total
}

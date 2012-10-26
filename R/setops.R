variantKeys <- function(x) paste0(x$location, x$alt)

setGeneric("%variant_in%", function(x, y) standardGeneric("%variant_in%"))

setMethod("%variant_in%", c("GenomicRanges", "GenomicRanges"),
          function(x, y) {
            hits <- findOverlaps(x, y)
            same.alt <- x$alt[queryHits(hits)] == y$alt[subjectHits(hits)]
            same <- logical(length(x))
            same[queryHits(hits)[same.alt]] <- TRUE
            same
          })

variant_setdiff <- function(x, y) {
  x[!x %variant_in% y]
}

## note this is asymmetric in that it keeps the metadata from 'x'
variant_intersect <- function(x, y) {
  x[x %variant_in% y]
}

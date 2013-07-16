
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Variant GRanges Compatibility
###

variantGRanges2VRanges <- function(x) {
  builtin.cols <- c("ref", "alt", "high.quality.ref", "high.quality",
                    "high.quality.total")
  mcols <- mcols(x)[setdiff(colnames(mcols(x)), builtin.cols)]
  with(mcols(x), VRanges(seqnames(x), ranges(x), ref, alt, high.quality.total,
                         high.quality.ref, high.quality, mcols))
}

variantGRangesIsDeprecated <- function(old) {
  .Deprecated("the VRanges class, see help(variantGRanges2VRanges)",
              old = old)
}


setMethod("ref", "GenomicRanges", function(x) x$ref)
setMethod("alt", "GenomicRanges", function(x) {
  if (is.null(x$alt))
    x$read
  else x$alt
})
setMethod("totalDepth", "GenomicRanges", function(x) {
  ifelse(is.na(x$high.quality.total), rawTotalDepth(x), x$high.quality.total)
})
setMethod("refDepth", "GenomicRanges", function(x) {
  ifelse(is.na(x$high.quality.ref), rawRefDepth(x), x$high.quality.ref)
})
setMethod("altDepth", "GenomicRanges", function(x) {
  ifelse(is.na(x$high.quality), rawAltDepth(x), x$high.quality)
})

## some internal compatibility wrappers for older versions
rawAltDepth <- function(x) {
  if (!is.null(x$raw.count))
    x$raw.count
  else x$count
}
rawTotalDepth <- function(x) {
  if (!is.null(x$raw.count.total))
    x$raw.count.total
  else x$count.total
}
rawRefDepth <- function(x) {
  if (!is.null(x$raw.count.ref))
    x$raw.count.ref
  else x$count.ref
}

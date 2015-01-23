### FIXME: this is too specific to the LR test. To make it more
### generic, we could combine p.error and p.lower into a special
### params object that would be used both for calling and for this
### function. Then, we would have a method that dispatches on that
### class. This params object could be a FilterRule, but we are
### currently not using polymorphism for FilterRule objects.

minCallableCoverage <- function(calling.filters, power = 0.80,
                                max.coverage = 1000L)
{
  if (!isSingleNumber(power) || power < 0 || power > 1)
    stop("'power' must be a single, non-NA number in [0,1]")
  if (!isSingleNumber(max.coverage) || max.coverage < 0)
    stop("'max.coverage' must be a single, non-negative, non-NA number")
  if (!is(calling.filters, "FilterRules"))
    stop("'calling.filters' must be a FilterRules object")
  lr.filter <- calling.filters$likelihoodRatio
  if (is.null(lr.filter))
    stop("'likelihoodRatio' filter not found in 'calling.filters'")
  rc.filter <- calling.filters$readCount
  if (is.null(rc.filter))
    stop("'readCount' filter not found in 'calling.filters'")
  size <- seq(1L, max.coverage)
  f <- lrtFreqCutoff(params(lr.filter)$p.lower, params(lr.filter)$p.error)
  p <- 1 - pbinom(round(pmax(params(rc.filter)$min.depth, size * f)), size,
                  params(lr.filter)$p.lower)
  cov <- head(size[p > power], 1L)
  if (length(cov) == 0L)
    NA_integer_
  else cov
}

setGeneric("callCallable", function(x, ...) standardGeneric("callCallable"))

setMethod("callCallable", "BigWigFile", function(x, pos = NULL, ...) {
  which <- if (is.null(pos)) seqinfo(x) else pos
  rle <- import(x, which = which, as="RleList")
  callCallable(rle, pos = pos, ...)
})

setMethod("callCallable", "RleList", function(x, pos = NULL, ...) {
  if (!is.null(pos)) {
    if (!is(pos, "GenomicRanges"))
      stop("'pos' must be NULL or a GenomicRanges")
    pos.all.width.one <- all(width(pos) == 1L)
    if (!pos.all.width.one)
      stop("All 'pos' ranges must be of length one")
  }
  cutoff <- minCallableCoverage(...)
  if (!is.null(pos)) {
    cov <- extractCoverageForPositions(x, pos)
    cov >= cutoff
  } else {
    x >= cutoff
  }
})

setMethod("callCallable", "ANY", function(x, ...) {
  cov <- coverage(x, drop.D.ranges = TRUE)
  callCallable(cov, ...)
})

# FIXME: instead of 'pos' have a 'which' and a logical byPos argument
# That way, people can iterate over partitions
callWildtype <- function(reads, variants, calling.filters, pos = NULL,
                         ...)
{
  callable <- callCallable(reads, calling.filters = calling.filters,
                           pos = pos, ...)
  wildtype <- callable
  callable[is.na(callable)] <- FALSE
  wildtype[!callable] <- NA
  if (is.null(pos)) {
    seqlevels(variants) <- names(callable)
    var <- as(variants, "RangesList")
  } else {
    var <- pos %over% variants
  }
  wildtype[var] <- FALSE
  wildtype
}

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
  f <- freq(params(lr.filter)$p.lower, params(lr.filter)$p.error, 1L)
  p <- 1 - pbinom(round(pmax(params(rc.filter)$min.count, size * f)), size,
                  params(lr.filter)$p.lower)
  head(size[p > power], 1L)
}

setGeneric("callCallable", function(x, ...) standardGeneric("callCallable"))

setMethod("callCallable", "RleList", function(x, param, ...) {
  if (!is(param, "BamTallyParam"))
    stop("'param' must be a BamTallyParam")
  cutoff <- minCallableCoverage(...)
  which.rl <- bamWhich(param)
  callable <- x < 0L # easy way to make FALSE RleList
  callable[which.rl] <- x[which.rl] >= cutoff
  callable
})
setMethod("callCallable", "ANY", function(x, ...) {
  cov <- coverage(x, drop.D.ranges = TRUE)
  callCallable(cov, ...)
})

callWildtype <- function(reads, variants, calling.filters, param, ...) {
  callable <- callCallable(reads, calling.filters, param, ...)
  wildtype <- callable
  wildtype[!callable] <- NA
  seqlevels(variants) <- names(callable)
  var.rl <- as(variants, "RangesList")
  wildtype[var.rl] <- FALSE
  wildtype
}

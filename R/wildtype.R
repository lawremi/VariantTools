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
  f <- freq(params(lr.filter)$p.lower, params(lr.filter)$p.error)
  p <- 1 - pbinom(round(pmax(params(rc.filter)$min.count, size * f)), size,
                  params(lr.filter)$p.lower)
  head(size[p > power], 1L)
}

setGeneric("callCallable", function(x, ...) standardGeneric("callCallable"))

setMethod("callCallable", "RleList", function(x, param, global = TRUE, ...) {
  if (!is(param, "BamTallyParam"))
    stop("'param' must be a BamTallyParam")
  cutoff <- minCallableCoverage(...)
  which.rl <- bamWhich(param)
  if (length(which.rl) > 0L) {
    if (global) {
      if (length(setdiff(seqlevels(which.rl), names(x))) > 0L)
        stop("Some seqlevels are missing from coverage")
      common.seqnames <- intersect(names(which.rl), names(x))
      x.in.which <- x[common.seqnames]
      callable <- seqapply(elementLengths(x.in.which), Rle, values = NA)
      callable[which.rl] <- x.in.which[which.rl] >= cutoff
      callable
    } else {
      cov <- extractCoverageForPositions(x, as(which.rl, "GRanges"))
      cov >= cutoff
    }
  } else {
    x >= cutoff
  }
})

setMethod("callCallable", "ANY", function(x, param, ...) {
  scan.param <- ScanBamParam(which = bamWhich(param))
  cov <- coverage(x, drop.D.ranges = TRUE, param = scan.param)
  callCallable(cov, param, ...)
})

callWildtype <- function(reads, variants, calling.filters, param, global = TRUE,
                         ...)
{
  if (!isTRUEorFALSE(global))
    stop("'global' must be TRUE or FALSE")
  which.all.width.one <-
    all(width(unlist(bamWhich(param), use.names = FALSE)) == 1L)
  if (!global && !which.all.width.one)
    stop("All 'bamWhich' ranges must be of length one for 'global = FALSE'")
  callable <- callCallable(reads, calling.filters = calling.filters,
                           param = param, global = global, ...)
  wildtype <- callable
  callable[is.na(callable)] <- FALSE
  wildtype[!callable] <- NA
  if (global) {
    seqlevels(variants) <- names(callable)
    var <- as(variants, "RangesList")
  } else {
    var <- bamWhich(param) %in% variants
  }
  wildtype[var] <- FALSE
  wildtype
}

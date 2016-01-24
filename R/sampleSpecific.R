### =========================================================================
### Finding Somatic Mutations
### -------------------------------------------------------------------------

setGeneric("callSampleSpecificVariants", function(case, control, ...) {
  standardGeneric("callSampleSpecificVariants")
})

### TODO: check for positions that are completely non-ref in normal
### and somehow gain the reference allele in tumor.

## NEW ALGORITHM
## - Filter out all normal calls from tumor calls
## - For the remaining tumor calls, ask two questions:
##   - Did we have sufficient coverage in normal to make a call?
##   - Is the normal frequency reasonable, given the tumor frequency?

## IDEAS FROM ROBERT:
## - Assume het frequency (~50%) for power and extremity filter.
##   - But this will miss artifacts in the normal, and thus
##     conclude that they are tumor-specific...
##   - The power filter might work best as an annotation.
## - Maybe even assume het frequency when calling in normal,
##   which would automatically fix the power filter.
## - Maybe discard tumor variants that look het, but risky;
##   analysts always have access to frequencies.

setMethod("callSampleSpecificVariants", c("VRanges", "VRanges"),
          function(case, control, control.cov, ...)
          {
            filters <- SampleSpecificVariantFilters(control,
                                                    control.cov,
                                                    hardFilters(case),
                                                    ...)
            case.specific <- subsetByFilter(case, filters)
            annotateWithControlDepth(case.specific, control, control.cov)
          })

setMethod("callSampleSpecificVariants", c("BamFile", "BamFile"),
          function(case, control, tally.param,
                   calling.filters = VariantCallingFilters(),
                   post.filters = FilterRules(),
                   ...)
          {
            case.raw <- tallyVariants(case, tally.param)
            case.called <- callVariants(case.raw, calling.filters, post.filters)

            control.raw <- tallyVariants(control, tally.param)
            sbp <- ScanBamParam(which = tally.param@bamTallyParam@which)
            control.cov <- coverage(control, drop.D.ranges = TRUE, param = sbp)

            callSampleSpecificVariants(case.called, control.raw,
                                       control.cov=control.cov,
                                       ...)
          })

setMethod("callSampleSpecificVariants", c("character", "character"),
          function(case, control, ...)
          {
            callSampleSpecificVariants(BamFile(case), BamFile(control), ...)
          })

SampleSpecificVariantFilters <-
  function(control, control.cov, calling.filters, power = 0.8, p.value = 0.01)
{
  control.called <- callVariants(control, calling.filters,
                                 post.filters = FilterRules())
  FilterRules(c(calledInControl = SetdiffVariantsFilter(control.called),
                power = CallableInOtherFilter(control.cov, calling.filters,
                  power),
                extremity = LowerFrequencyInOtherFilter(control, control.cov,
                  p.value)
                ))
}

SetdiffVariantsFilter <- function(other) {
  function(x) {
    !(x %in% other)
  }
}

SetdiffPositionsFilter <- function(other) {
  function(x) {
    !(x %over% other)
  }
}

extractCoverageForPositions <- function(cov, pos) {
  if (length(setdiff(seqlevels(pos), names(cov))) > 0L)
    stop("Some seqlevels are missing from coverage")
  if (any(width(pos) > 1L))
    stop("Some ranges are of width > 1")
  seqlevels(pos) <- names(cov)
  ord <- order(seqnames(pos))
  ans <- integer(length(pos))
  ans[ord] <- unlist(mapply(function(v, p) {
    runValue(v)[findRun(p, v)]
  }, cov, split(start(pos), seqnames(pos)), SIMPLIFY=FALSE), use.names=FALSE)
  ans
}

CallableInOtherFilter <-
  function(other.cov, calling.filters, min.power = 0.8,
           p.lower = params(calling.filters$likelihoodRatio)$p.lower)
{
  function(x) {
    calculatePowerInOther(x, other.cov, calling.filters, p.lower) >= min.power
  }
}

calculatePowerInOther <-
  function(x, other.cov, calling.filters = hardFilters(x),
           p.lower = params(calling.filters$likelihoodRatio)$p.lower)
{
  lr.filter <- calling.filters$likelihoodRatio
  if (is.null(lr.filter))
    f <- 0
  else
    f <- lrtFreqCutoff(params(lr.filter)$p.error, params(lr.filter)$p.lower)
  rc.filter <- calling.filters$readCount
  if (is.null(rc.filter))
    min.depth <- 0L
  else
    min.depth <- params(rc.filter)$min.depth
  other.n <- extractCoverageForPositions(other.cov, resize(x, 1))
  min.depth <- ceiling(pmax(min.depth, other.n * f))
  1 - pbinom(min.depth-1L, other.n, p.lower)
}

LowerFrequencyInOtherFilter <- function(other, other.cov, p.value = 0.01)
{
  function(x) {
      ## NOTE: this was once based on the raw (not quality filtered)
      ## counts, but bam_tally is no longer giving us those for the alts.
    x <- annotateWithControlDepth(x, other, other.cov)
    p <- with(x, pbinom(control.alt.depth, control.total.depth,
                        altDepth / totalDepth))
    p < p.value
  }
}

annotateWithControlDepth <- function(case, control, control.cov) {
  m <- match(case, control)
  control.alt.depth <- as.vector(altDepth(control))[m]
  control.alt.depth[is.na(control.alt.depth)] <- 0L
  control.total.depth <- as.vector(totalDepth(control))[m]
  control.raw.total.depth <- as.vector(rawTotalDepth(control))[m]
  control.raw.total.depth[is.na(m)] <-
      extractCoverageForPositions(control.cov, resize(case[is.na(m)], 1))
  control.alt.depth[is.na(m)] <- control.raw.total.depth[is.na(m)]
  case$control.alt.depth <- control.alt.depth
  case$control.total.depth <- control.total.depth
  case$control.raw.total.depth <- control.raw.total.depth
  case
}

caseControlFET <- function(case, control, control.cov) {
  case <- annotateWithControlDepth(case, control, control.cov)
  with(case,
       fisher_p(altDepth, (totalDepth - altDepth),
                altDepth + control.alt.depth,
                (totalDepth - altDepth) +
                (control.total.depth - control.alt.depth)))
}

DepthFETFilter <- function(control, control.cov, p.value.cutoff = 0.05) {
  function(x) {
    p.value <- caseControlFET(x, control, control.cov)
    p.value < p.value.cutoff
  }
}

MaxControlFreqFilter <- function(control, control.cov, max.control.freq = 0.03) {
  function(x) {
    x <- annotateWithControlDepth(x, control, control.cov)
    freq <-  x$control.alt.depth / x$control.total.depth
    ifelse(is.na(freq), 0, freq) <= max.control.freq
  }
}

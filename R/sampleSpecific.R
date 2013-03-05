### =========================================================================
### Finding Somatic Mutations
### -------------------------------------------------------------------------

## For the API, we need to decide how this fits into the workflow.
## For example, we could go directly from the BAMs, tallying and
## calling the variants, etc. Or we could start from the
## sample-specific calls and the normal raw tallies. And eventually we
## will also need to get at the tumor raw reference tallies, probably
## dipping back into the BAM. So the API is more complicated in the
## second case, but much time is saved by not needing to repeat the
## tally. There is another issue: for the power calculation we are
## assuming that the calls were generated according to the same method
## and parameters. In theory, the method, parameters, BAM, etc, could
## be attached to the objects in their metadata. There are no
## guarantees, but it would be a convenient feature. Another option is
## to take the raw variants for both tumor and normal, and perform the
## sample-specific calling first. This is way less error-prone than
## the previous approach. That would be a little slow, due to the
## fisher test, which is always slow (do we need some C?). What about
## the coverage though? Do we allow the user to pass it? Is it
## calculated from a BAM in the metadata? For now, user passes it.

## Anyway,
## this is just the high-level API. The low-level API is always
## available for people that want to play around.

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

## Argument for no QA of control: if a variant is filtered out by QA,
## and it passes the calling filters, we do not feel comfortable
## saying that the variant does not exist in control. Most of the
## time, we will recover the variant call by the extremity
## test. However, if the case variant is present at a high fraction
## (possbibly due to LOH), and the control variant is still callable,
## but much lower in fraction, we may miss it. One way to get around
## that would be to use 'plower' as the cap for the case
## frequency. But why throw them out by QA in the first place? This
## really comes down to a lack of confidence in our QA filters: they
## are meant to flag/remove questionable data; not confirm that there
## really is a problem.

## Question: in the extremity test, do we have to consider zeros ALT
## counts in control? We will not have those in the tallies. Does the
## power check cover it? Coverage may indeed be high enough to make
## the call, while the extremity test could still fail. Maybe.

setMethod("callSampleSpecificVariants", c("GenomicRanges", "GenomicRanges"),
          function(case, control, control.cov,
                   qa.filters = VariantQAFilters(),
                   calling.filters = VariantCallingFilters(),
                   post.filters = VariantPostFilters(),
                   ...)
          {
            case.called <- callVariants(qaVariants(case, qa.filters),
                                        calling.filters, post.filters)
            
            filters <- SampleSpecificVariantFilters(control,
                                                    control.cov,
                                                    calling.filters,
                                                    ...)
            case.specific <- subsetByFilter(case.called, filters)
            annotateWithControlCounts(case.specific, control, control.cov)
          })

setMethod("callSampleSpecificVariants", c("BamFile", "BamFile"),
          function(case, control, tally.param, ...)
          {
            case.raw <- tallyVariants(case, tally.param)
            control.raw <- tallyVariants(control, tally.param)

            control.cov <- coverage(control, drop.D.ranges = TRUE)
            
            callSampleSpecificVariants(case.raw, control.raw, control.cov, ...)
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
    !(x %variant_in% other)
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
  ord <- order(pos)
  pos <- pos[ord]
  rl <- as(pos, "RangesList")
  ans <- integer(length(pos))
  ans[ord] <- as.vector(unlist(seqselect(cov, rl), use.names = FALSE))
  ans
}

CallableInOtherFilter <- function(other.cov, calling.filters, power = 0.8)
{
  lr.filter <- calling.filters$likelihoodRatio
  if (is.null(lr.filter))
    stop("'likelihoodRatio' filter not found in 'calling.filters'")
  lr.params <- params(lr.filter)
  rc.filter <- calling.filters$readCount
  if (is.null(rc.filter))
    stop("'readCount' filter not found in 'calling.filters'")
  function(x) {
    other.n <- extractCoverageForPositions(other.cov, x)
    f <- lrtFreqCutoff(lr.params$p.error, lr.params$p.lower)
    min.count <- round(pmax(params(rc.filter)$min.count, other.n * f))
    1 - pbinom(min.count, other.n, lr.params$p.lower) > power
  }
}

LowerFrequencyInOtherFilter <- function(other, other.cov, p.value = 0.01)
{
  function(x) {
    x.freq <- with(mcols(x), high.quality / high.quality.total)
    m <- match(variantKeys(x), variantKeys(other))
    other.alt <- rep.int(0L, length(x))
    other.alt[!is.na(m)] <- other$count[m[!is.na(m)]]
    other.total <- extractCoverageForPositions(other.cov, x)
    p <- pbinom(other.alt, other.total, x.freq)
    p < p.value
  }
}

annotateWithControlCounts <- function(case.specific, control, control.cov) {
  m <- matchVariants(case.specific, control)
  control.count <- control$high.quality[m]
  control.count[is.na(control.count)] <- 0L
  control.count.total <- control$high.quality.total[m]
  control.count.total[is.na(m)] <-
    extractCoverageForPositions(control.cov, case.specific[is.na(m)])
  case.specific$control.count <- control.count
  case.specific$control.count.total <- control.count.total
  case.specific
}

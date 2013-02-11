### =========================================================================
### Variant Calling
### -------------------------------------------------------------------------

setGeneric("callVariants", function(x, ...) standardGeneric("callVariants"))

## This method has explicit arguments for each stage in the
## pipeline. This is a concious choice to expose the algorithmic
## details (but not the implementation details).

setMethod("callVariants", "BamFile",
          function(x, tally.param,
                   qa.filters = VariantQAFilters(),
                   calling.filters = VariantCallingFilters(...),
                   post.filters = VariantPostFilters(),
                   ...)
          {
            raw_variants <- tallyVariants(x, tally.param)
            qa_variants <- qaVariants(raw_variants, qa.filters)
            callVariants(qa_variants, calling.filters, post.filters)
          })

setMethod("callVariants", "character", function(x, ...) {
  callVariants(BamFile(x), ...)
})

setMethod("callVariants", "GenomicRanges",
          function(x, calling.filters = VariantCallingFilters(...),
                   post.filters = VariantPostFilters(), ...)
          {
            called_variants <- subsetByFilter(x, calling.filters)
            postFilterVariants(called_variants, post.filters)
          })

VariantCallingFilters <-
  function(read.count = 2L, p.lower = 0.2, p.error = 1/10000)
{
  c(VariantSanityFilters(),
    FilterRules(list(readCount = MinCountFilter(read.count),
                     likelihoodRatio = BinomialLRFilter(p.lower, p.error))))
}

MinCountFilter <- function(min.count = 2L) {
  function(x) {
    mcols(x)[["high.quality"]] >= min.count
  }
}

BinomialLRFilter <-
  function(p.lower = 0.2, p.error = 1/10000)
{
  function(x) {
    freq <- freq(p.error, p.lower)
    sampleFreq <- with(values(x), high.quality / (high.quality.total))
    passed <- sampleFreq >= freq
    passed[is.na(passed)] <- FALSE
    passed
  }
}

freq <- function(p0, p1, n = if (C == 1L) 1L else n, C = 1L) {
  num <- (1/n) * log(C) + log(1-p0) - log(1-p1)
  denom <- log(p1) - log(p0) + log(1-p0) - log(1-p1)
  num/denom
}

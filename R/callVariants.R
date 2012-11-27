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
                   ...)
          {
            raw_variants <- tallyVariants(x, tally.param)
            qa_variants <- qaVariants(raw_variants, qa.filters)
            callVariants(qa_variants, calling.filters)
          })

setMethod("callVariants", "character", function(x, ...) {
  callVariants(BamFile(x), ...)
})

setMethod("callVariants", "GenomicRanges",
          function(x, calling.filters = VariantCallingFilters(...), ...)
          {
            subsetByFilter(x, calling.filters)
          })

VariantCallingFilters <-
  function(read.count = 2L, p.lower = 0.2,
           p.error = if (use.high.qual) 1/10000 else 1/100,
           use.high.qual = TRUE)
{
  c(VariantSanityFilters(),
    FilterRules(list(readCount = MinCountFilter(read.count, use.high.qual),
                     likelihoodRatio =
                     BinomialLRFilter(p.lower, p.error, use.high.qual))))
}

MinCountFilter <- function(min.count = 2L, use.high.qual = TRUE) {
  function(x) {
    mcols(x)[[if (use.high.qual) "high.quality" else "count"]] >= min.count
  }
}

BinomialLRFilter <-
  function(p.lower = 0.2,
           p.error = if (use.high.qual) 1/10000 else 1/100,
           use.high.qual = TRUE)
{
  function(x) {
    freq <- freq(p.error, p.lower, 1)
    sampleFreq <- if (use.high.qual) {
      with(values(x), high.quality / (high.quality.total))
    } else {
      with(values(x), count / count.total)
    }
    passed <- sampleFreq >= freq
    passed[is.na(passed)] <- FALSE
    passed
  }
}

freq <- function(p0, p1, n, C) {
  num <- (1/n) * log(C) + log(1-p0) - log(1-p1)
  denom <- log(p1) - log(p0) + log(1-p0) - log(1-p1)
  num/denom
}

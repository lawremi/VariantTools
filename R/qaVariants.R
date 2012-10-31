### =========================================================================
### QA Variant Filtering (as opposed to the calling filters)
### -------------------------------------------------------------------------

sanitizeVariants <- function(x, ...) {
  subsetByFilter(x, VariantSanityFilters(...))
}

VariantSanityFilters <- function() {
  FilterRules(list(nonRef = NonRefFilter(), nonNRef = NonNRefFilter()))
}

qaVariants <- function(x, qa.filters = VariantQAFilters(...), ...)
{
  subsetByFilter(x, qa.filters)
}

VariantQAFilters <- function(cycle.count = 2L, fisher.strand.p.value = 1e-3)
{
  c(VariantSanityFilters(),
    FilterRules(c(cycleCount = CycleCountFilter(cycle.count),
                  fisherStrand = FisherStrandFilter(fisher.strand.p.value),
                  cycleBin = InternalCycleBinFilter())))
}

## With new gmapR, this is only necessary for filtering the ref N's.
## In theory, we could keep positions with ref N and at least two alts.
## But this is OK for now.
NonNRefFilter <- function() {
  function(x) {
    as.character(x$ref) != 'N'
  }
}

## Drops the ref rows (should not be necessary)
NonRefFilter <- function() {
  function(x) {
    !is.na(x$alt)
  }
}

CycleCountFilter <- function(cycle.count = 2L) {
  function(x) {
    x$ncycles >= cycle.count
  }
}

BinomialErrorFilter <- function(p.error = 1/1000, p.value = 0.01) {
  function(x) {
    p <- pbinom(values(x)$count - 1, values(x)$count.total, prob = p.error,
                lower.tail = FALSE)
    p < p.value
  }
}

FisherStrandFilter <- function(p.value = 1e-3) {
  function(x) {
    p <- with(mcols(x),
              fisher_p_vectorized(count.pos.ref,
                                  count.neg.ref,
                                  (count.pos.ref + count.pos),
                                  (count.neg.ref + count.neg)))
    p > p.value
  }
}

InternalCycleBinFilter <- function(min.count = 1L) {
  function(x) {
    cycle_columns <- grep("cycleCount", colnames(mcols(x)), value = TRUE)
    alt_columns <- grep("ref", cycle_columns, invert = TRUE, value = TRUE)
    internal_columns <- tail(head(alt_columns, -1), -1)
    if (length(internal_columns) > 0L)
      rowSums(as.matrix(mcols(x)[internal_columns])) >= min.count
    else rep.int(TRUE, length(x))
  }
}

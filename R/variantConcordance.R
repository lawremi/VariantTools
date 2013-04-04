
calculateVariantConcordance <- function(gr1, gr2, which = NULL) {
  filters <- VariantConcordanceFilters()
  
  gr1 <- subsetByFilter(gr1, filters)
  gr2 <- subsetByFilter(gr2, filters)
  
  if(!is.null(which)){
    gr1 <- subsetByOverlaps(gr1, which)
    gr2 <- subsetByOverlaps(gr2, which)
  }
  
  overlaps <- findOverlaps(gr1, gr2)
  test1 <-gr1[queryHits(overlaps)]
  test2 <-gr2[subjectHits(overlaps)]

  callHetHom <- function(x) {
    alt <- if (is.null(x$alt)) x$read else x$alt
    ifelse(x$count.ref > 2,
           paste0(as.character(x$ref), as.character(alt)),
           paste0(as.character(alt), as.character(alt)))
  }

  t1geno <- callHetHom(test1)
  t2geno <- callHetHom(test2)
  
  fraction <- sum(t1geno == t2geno) / length(test1)

  DataFrame(fraction = fraction, denominator = length(test1))
}

SingleVariantFilter <- function() {
  function(x) {
    location <- Rle(x$location)
    rep(width(location) == 1L, width(location))
  }
}

VariantConcordanceFilters <- function() {
  FilterRules(c(NonN = NonNRefFilter(),
                minCount = MinCountFilter(10L),
                singleVariant = SingleVariantFilter()))
}

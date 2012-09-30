
checkVariantConcordance <- function(gr1, gr2, SNP_gr = NULL) {
  filters <- VariantConcordanceFilters()
  
  gr1 <- subsetByFilter(gr1, filters)
  gr2 <- subsetByFilter(gr2, filters)
  
  if(!is.null(SNP_gr)){
    gr1 <- subsetByOverlaps(gr1, SNP_gr)
    gr2 <- subsetByOverlaps(gr2, SNP_gr)
  }
  
  overlaps <- findOverlaps(gr1, gr2)
  test1 <-gr1[queryHits(overlaps)]
  test2 <-gr2[subjectHits(overlaps)]

  callHetHom <- function(x) {
    ifelse(x$count.ref > 2,
           paste0(as.character(x$ref), as.character(x$alt)),
           paste0(as.character(x$alt), as.character(x$alt)))
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

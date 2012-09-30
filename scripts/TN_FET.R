TwoSampleFET <- function(gr, raw, cov, pval=0.001, freq=0.05){
  values(gr)$norm_ref <- as.numeric(NA)
  values(gr)$norm_alt <- as.numeric(NA)
  over2 <- which(paste(values(gr)$location,values(gr)$read) %in% paste(values(raw)$location, values(raw)$read)) 
  over <- which(paste(values(raw)$location,values(raw)$read) %in% paste(values(gr)$location, values(gr)$read))
  values(gr[over2])$norm_ref <- values(raw[over])$count.ref
  values(gr[over2])$norm_alt <- values(raw[over])$count
  values(gr[-over2])$norm_alt <- 0

  grl <- split(gr, as.factor(seqnames(gr)))
  names <- as.list(names(grl))
  n <- lapply(names, function(x) { as.numeric(cov[[x]][ranges(grl[[x]])]) })
  n <- unlist(n)
  gr <- unlist(grl)
  values(gr)$norm_cov<-n
  pvals <- fisher_p(values(gr)$count, (values(gr)$count.total -values(gr)$count),
           values(gr)$count+values(gr)$norm_alt,
           (values(gr)$count.total-values(gr)$count)+(values(gr)$norm_cov-values(gr)$norm_alt))
  values(gr)$pval <- pvals
  TS <- gr[values(gr)$pval<pval &
                 values(gr)$norm_alt/
                 (values(gr)$norm_cov)<freq]
  return(TS)
}

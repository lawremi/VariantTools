callableFraction <- function(tx, reads, ...) {
  callable <- callCallable(reads, ...)
  seqlevels(tx) <- names(cov)
  tx_gr <- unlist(tx, use.names = FALSE)
  tx_ranges <- split(ranges(tx_gr),
                     seqnames(tx_gr), names(cov))
  rowsum(unlist(viewSums(Views(tx_ranges, callable)), use.names = FALSE),
         names(tx)[togroup(tx)]) / sum(width(tx))[togroup(tx)]
}

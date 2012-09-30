
### FIXME: this is too specific to the LR test. To make it more
### generic, we could combine p.error and p.lower into a special
### params object that would be used both for calling and for this
### function. Then, we would have a method that dispatches on that
### class. This params object could be a FilterRule, but we are
### currently not using polymorphism for FilterRule objects.

callCallable <- function(reads, plower = 0.2, perror = 1/1000, power = 0.999)
{
  ## would it be cleaner to find all unique coverage values and
  ## use that for 'size'?
  size <- seq(20, 1000)
  f <- freq(n = 1, p1 = plower, p0 = perror)
  p <- 1 - pbinom(round(size * f), size, plower)
  ## what happens when which() is empty?
  cutoff <- min(which(p > power)) + 20
  cov <- coverage(reads, drop.D.ranges = TRUE)
  cov > cutoff
}

callableFraction <- function(tx, reads, ...) {
  callable <- callCallable(reads, ...)
  seqlevels(tx) <- names(cov)
  tx_gr <- unlist(tx, use.names = FALSE)
  tx_ranges <- split(ranges(tx_gr),
                     seqnames(tx_gr), names(cov))
  rowsum(unlist(viewSums(Views(tx_ranges, callable)), use.names = FALSE),
         names(tx)[togroup(tx)]) / sum(width(tx))[togroup(tx)]
}

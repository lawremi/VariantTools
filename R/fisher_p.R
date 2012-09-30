fisher_p <- function( x, y, m, n, relErr = 1 + 1e-7 ) {
   # 'x' and 'y' are entries in the top two cells; 'm' and 'n' are
   # column totals.  Code is excerpted from fisher.test, for
   # efficiency. Note that 'support' varies in length with the input
   # variables, so vectorization is only possible via an mapply.
   mapply(
          function( x, y, m, n ) {
            k <- x + y
            lo <- max( 0, k - n )
            hi <- min( k, m )
            support <- lo:hi
            d <- dhyper( support, m, n, k, log = TRUE )
            d <- exp( d - max( d ) )
            d <- d / sum( d )
            sum( d[ d <= d[ x - lo + 1 ] * relErr ] )
          },
          x, y, m, n
          )
 }

fisher_p_vectorized <- function(x, y, m, n, relErr = 1 + 1e-7) {
  ## 'x' and 'y' are entries in the top two cells; 'm' and 'n' are
  ## column totals. Vectorized version of Richard's fisher_p().
  ## Only about 3X faster than non-vectorized version.
  ## With some simple sum and max operations over partitions,
  ## this would be way faster.
  k <- x + y
  lo <- pmax(0L, k - n)
  hi <- pmin(k, m)
  support <- IRanges:::mseq(lo, hi)
  w <- hi - lo + 1L
  m <- rep(m, w)
  n <- rep(n, w)
  k <- rep(k, w)
  d <- dhyper(support, m, n, k, log = TRUE)
  part <- PartitioningByWidth(w)
  group <- togroup(part)
  maxd <- viewMaxs(Views(Rle(d), part))
  ##maxd <- d[order(group, d)[end(part)]]
  d <- exp(d - rep(maxd, w))
  d <- d / rep(rowsum(d, group, reorder = FALSE), w)
  sel <- d <= rep(d[x - lo + start(part)], w) * relErr
  unname(rowsum(d[sel], group[sel], reorder = FALSE)[,1])
}

## Unused attempt using IRanges Lists. Not much of a speed improvement.
fisher_p_List <- function( x, y, m, n, relErr = 1 + 1e-7 ) {
  k <- x + y
  lo <- pmax(0L, k - n)
  hi <- pmin(k, m)
  ranges <- IRanges(lo, hi)
  support <- as.integer(ranges)
  m <- rep(m, width(ranges))
  n <- rep(n, width(ranges))
  k <- rep(k, width(ranges))
  d <- relist(Rle(dhyper(support, m, n, k, log = TRUE)),
              PartitioningByWidth(width(ranges)))
  d <- exp(d - as(max(d), "List"))
  d <- d / as(sum(d), "List")
  sum(d[d <= d[as(x - lo + 1L, "List")] * relErr])
}

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
  k <- x + y
  lo <- pmax(0L, k - n)
  hi <- pmin(k, m)
  support <- IRanges:::unlist_as_integer(IRanges(lo, hi))
  w <- hi - lo + 1L
  m <- rep(m, w)
  n <- rep(n, w)
  k <- rep(k, w)
  d <- dhyper(support, m, n, k, log = TRUE)
  part <- PartitioningByWidth(w)
  dl <- relist(d, part)  
  maxd <- max(dl)
  dl <- exp(dl - maxd) / sum(dl)
  sel <- dl <= dl[as(x - lo + 1L, "IntegerList")] * relErr
  unname(sum(dl[sel]))
}

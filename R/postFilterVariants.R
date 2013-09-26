### =========================================================================
### Post filtering of variants (come after the calling filters)
### -------------------------------------------------------------------------

postFilterVariants <- function(x, post.filters = VariantPostFilters(...), ...) {
  subsetByFilter(x, post.filters)
}

VariantPostFilters <- function(max.nbor.count = 0.1, whitelist = NULL) {
  FilterRules(list(avgNborCount = AverageNeighborCountFilter(max.nbor.count,
                     whitelist = whitelist)))
}

twoWayTabulate <- function(x, y, x.nbins, y.nbins) {
  in.bin <- x > 0L & x <= x.nbins & y > 0L & y <= y.nbins
  new("dgTMatrix", i = x[in.bin] - 1L, j = y[in.bin] - 1L,
      x = rep.int(1, sum(in.bin)), Dim = c(x.nbins, y.nbins))
}

findAverageNeighborCount <- function(x, max.dist = 50L,
                                     weighter = function(d) 1/sqrt(d),
                                     neighbors = x)
{
  windows <- suppressWarnings(resize(neighbors, max.dist * 2L + 1L,
                                     fix = "center"))
  hits <- findOverlaps(x, windows)
  dist <- abs(start(neighbors)[subjectHits(hits)] - start(x)[queryHits(hits)])
  dist.counts <- twoWayTabulate(queryHits(hits), dist,
                                queryLength(hits), max.dist)
  weights <- weighter(seq_len(ncol(dist.counts)))
  norm.weights <- weights / sum(weights)
  suppressMessages(dist.counts %*% matrix(norm.weights))[,1]
}

### FIXME: do we want to consider the alt here, i.e.,
### should whitelist be a VRanges?
AverageNeighborCountFilter <- function(max.nbor.count = 0.1, max.dist = 50L,
                                       weighter = function(d) 1/sqrt(d),
                                       whitelist = NULL)
{
  function(x) {
    if (!is.null(whitelist)) {
      white.listed <- !MaskFilter(whitelist)(x)
      filtered.x <- x[!white.listed]
    } else filtered.x <- x
    count <- findAverageNeighborCount(filtered.x, max.dist, weighter)
    ans <- count <= max.nbor.count
    if (!is.null(whitelist)) {
      ans.padded <- rep.int(TRUE, length(x))
      ans.padded[!white.listed] <- ans
      ans.padded
    } else ans
  }
}

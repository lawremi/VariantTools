### =========================================================================
### Post filtering of variants (come after the calling filters)
### -------------------------------------------------------------------------

postFilterVariants <- function(x, post.filters = VariantPostFilters(...), ...) {
  subsetByFilter(x, post.filters)
}

VariantPostFilters <- function(nbor.max.count = 0.1) {
  FilterRules(list(avgNborCount = AverageNeighborCountFilter(nbor.max.count)))
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

AverageNeighborCountFilter <- function(nbor.max.count = 0.1, max.dist = 50L,
                                       weighter = function(d) 1/sqrt(d))
{
  function(x) {
    count <- findAverageNeighborCount(x, max.dist, weighter)
    count <= nbor.max.count
  }
}

### =========================================================================
### Variant QA (as opposed to the calling filters)
### -------------------------------------------------------------------------

qaVariants <- function(x, qa.filters = VariantQAFilters(...), ...)
{
  softFilter(x, qa.filters)
}

VariantQAFilters <- function(fisher.strand.p.value = 1e-4, min.mdfne = 10L)
{
  FilterRules(c(mdfne = MedianDistFromNearestEndFilter(min.mdfne),
                fisherStrand = FisherStrandFilter(fisher.strand.p.value)))
}

## With new gmapR, this is only necessary for filtering the ref N's.
## In theory, we could keep positions with ref N and at least two alts.
## But this is OK for now.
NonNRefFilter <- function() {
  function(x) {
    as.character(ref(x)) != 'N'
  }
}

## Drops the ref rows (should not be necessary)
NonRefFilter <- function() {
  function(x) {
    !is.na(alt(x))
  }
}

ReadPosCountFilter <- function(read.pos.count = 1L) {
  function(x) {
    (if (!is.null(x$ncycles)) x$ncycles else x$n.read.pos) >= read.pos.count
  }
}

FisherStrandFilter <- function(p.value = 1e-4) {
  function(x) {
    p <- with(mcols(x),
              fisher_p_vectorized(count.plus.ref,
                                  count.minus.ref,
                                  (count.plus.ref + count.plus),
                                  (count.minus.ref + count.minus)))
    p > p.value
  }
}

InternalReadPosBinFilter <- function(min.count = 1L) {
  function(x) {
    read_pos_columns <- grep("readPosCount", colnames(mcols(x)), value = TRUE)
    alt_columns <- grep("ref", read_pos_columns, invert = TRUE, value = TRUE)
    internal_columns <- tail(head(alt_columns, -1), -1)
    if (length(internal_columns) > 0L)
      rowSums(as.matrix(mcols(x)[internal_columns])) >= min.count
    else rep.int(TRUE, length(x))
  }
}

t.test_welch <- function(m1, m2, s1, s2, n1, n2) {
  s <- s1 / n1 + s2 / n2
  t <- (m1 - m2) / sqrt(s)
  v <- s^2 / ((s1 / n1)^2 / (n1 - 1L) + (s2 / n2)^2 / (n2 - 1L))
  pt(-abs(t), v) * 2L
}

ReadPositionTTestFilter <- function(p.value.cutoff = 1e-4) {
  function(x) {
    p <- with(mcols(x), t.test_welch(read.pos.mean, read.pos.mean.ref,
                                     read.pos.var, read.pos.var.ref,
                                     rawAltDepth(x), rawTotalDepth(x)))
    ans <- p > p.value.cutoff
    ans[is.na(ans)] <- TRUE
    ans
  }
}

DistanceToNearestFilter <- function(min.dist = 10L) {
  function(x) {
    mcols(distanceToNearest(x))$distance > min.dist
  }
}

NeighborCountFilter <- function(max.count = 2L, window.size = 75L) {
  function(x) {
    countOverlaps(resize(x, window.size, fix = "center"), x) <= max.count
  }
}

IndelsNotSupportedFilter <- function() {
  function(x) {
    nzchar(ref(x)) & nzchar(alt(x))
  }
}

MaskFilter <- function(mask) {
  function(x) {
    !overlapsAny(x, mask, ignore.strand = TRUE)
  }
}

MedianDistFromNearestEndFilter <- function(min.mdfne) {
  function(x) {
    x$mdfne >= min.mdfne
  }
}

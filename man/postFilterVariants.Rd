\name{postFilterVariants}
\alias{postFilterVariants}
\alias{VariantPostFilters}

\title{
  Post-filtering of Variants
}
\description{
  Applies filters to a set of called variants. The only current filter
  is a cutoff on the weighted neighbor count of each variant. This
  filtering is performed automatically by \code{\link{callVariants}}, so
  these functions are for when more control is desired.
}
\usage{
postFilterVariants(x, post.filters = VariantPostFilters(...), ...)
VariantPostFilters(max.nbor.count = 0.1, whitelist = NULL)
}

\arguments{
  \item{x}{
    A tally \code{GRanges} containing called variants, as output by
    \code{\link{callVariants}}.
  }
  \item{post.filters}{
    The filters applied to the called variants.
  }
  \item{\dots}{
    Arguments passed to \code{VariantPostFilters}, listed below.
  }
  \item{max.nbor.count}{
    Maximum allowed number of neighbors (weighted by distance)
  }
  \item{whitelist}{
    Positions to ignore; these will always pass the filter, and are
    excluded from the neighbor counting.
  }
}
\details{
  The neighbor count is calculated within a 100bp window centered on the
  variant. Each neighbor is weighted by the inverse square root of the
  distance to the neighbor. This was motivated by fitting logistic regression
  models including a term the count (usually 0, 1, 2) at each
  distance. The inverse square root function best matched the trend in
  the coefficients.
}
\value{
  For \code{postFilterVariants}, a tally \code{GRanges} of the variants that
  pass the filters.

  For \code{VariantPostFilters}, a \code{\link[IRanges]{FilterRules}}
  object with the filters.
}

\author{
  Michael Lawrence and Jeremiah Degenhardt
}

\examples{
p53 <- gmapR:::exonsOnTP53Genome("TP53")
bams <- LungCancerLines::LungCancerBamFiles()
tally.param <- VariantTallyParam(gmapR::TP53Genome(), 
                                 readlen = 100L,
                                 high_base_quality = 23L,
                                 which = range(p53))
# disable post-filtering during variant calling
called.variants <- callVariants(bam, tally.param,
                                post.filters = FilterRules())
# and apply at a later time...
postFilterVariants(called.variants, max.nbor.count = 0.15)
}

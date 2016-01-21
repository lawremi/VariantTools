\name{qaVariants}
\alias{qaVariants}
\alias{VariantQAFilters}

\title{
  QA Filtering of Variants
}
\description{
  Filters a tally \code{GRanges} through a series of simple checks for
  strand and read position (read position) biases.
}
\usage{
qaVariants(x, qa.filters = VariantQAFilters(...), ...)
VariantQAFilters(fisher.strand.p.value = 1e-4, min.mdfne = 10L)
}

\arguments{
  \item{x}{
    A tally \code{GRanges} as output by \code{\link{tallyVariants}}.
  }
  \item{qa.filters}{
    The filters used for the QA process, typically constructed with
    \code{VariantQAFilters}, see arguments below.
  }
  \item{\dots}{
    Arguments passed to \code{VariantQAFilters}, listed below.
  }
  \item{fisher.strand.p.value}{p-value cutoff for the Fisher's Exact Test for
    strand bias (+/- counts, alt vs. ref). Any variants with p-values
    below this cutoff are discarded.
  }
  \item{min.mdfne}{Minimum allowed median distance of alt calls from
    their nearest end of the read.
  }
}
\details{
  There are currently two QA filters:
  \itemize{
    \item{Median distance of alt calls from nearest end of the read is
      required to be >= \code{min.mdfne}, which defaults to 10.
    }
    \item{Fisher's Exact Test for strand bias, using the +/- counts, alt
      vs. ref. If the null is rejected, the variant is discarded.
    }
  }
}
\value{
  For \code{qaVariants}, a tally \code{GRanges} of the variants that
  pass the QA checks.

  For \code{VariantQAFilters}, a \code{\link[S4Vectors]{FilterRules}}
  object with the QA and sanity filters.
}

\author{
  Michael Lawrence and Jeremiah Degenhardt
}

\examples{
bams <- LungCancerLines::LungCancerBamFiles()
tally.param <- TallyVariantsParam(gmapR::TP53Genome(), 
                                  high_base_quality = 23L,
                                  which = gmapR::TP53Which())
tally.variants <- tallyVariants(bams$H1993, tally.param)
qaVariants(tally.variants, fisher.strand.p.value = 1e-4)
}

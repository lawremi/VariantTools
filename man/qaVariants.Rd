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
VariantQAFilters(read.pos.count = 1L, fisher.strand.p.value = 1e-4,
                 read.pos.p.value = 1e-4)
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
  \item{read.pos.count}{Minimum number of unique read positions for the
    alternate base.
  }
  \item{fisher.strand.p.value}{p-value cutoff for the Fisher's Exact Test for
    strand bias (+/- counts, alt vs. ref). Any variants with p-values
    below this cutoff are discarded.
  }
  \item{read.pos.p.value}{
    p-value cutoff for the read position t-test
    between the variant and reference calls
  }
}
\details{
  There are currently three QA filters:
  \itemize{
    \item{Alternate base was read at a minimum (1) number of unique
      read positions. We recommend setting this to 2 or more if
      duplicate alignments have not been filtered during tallying. 
    }
    \item{Fisher's Exact Test for strand bias, using the +/- counts, alt
      vs. ref. If the null is rejected, the variant is discarded.
    }
    \item{If the tallies contain read position bin counts, the variant must have
      at least one count in the middle bins
      (those not at the start or end). We trust the internal read positions more.
    }
    \item{Read position t-test comparing the mean read position for the
      reference and alt reads. Any imbalance probably indicates mapping
      issues. 
    }
    \item{Mask of blacklisted positions, such as simple repeats, low
      complexity regions, i.e., uninteresting, problematic regions.
    }
  }
  
  Prior to the QA checks, the variants are passed through a simple
  sanity filter that discards positions where reference has an N.
}
\value{
  For \code{qaVariants}, a tally \code{GRanges} of the variants that
  pass the QA checks.

  For \code{VariantQAFilters}, a \code{\link[IRanges]{FilterRules}}
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

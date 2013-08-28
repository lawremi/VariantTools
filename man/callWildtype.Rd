\name{callWildtype}
\alias{callWildtype}
\alias{minCallableCoverage}
\title{
  Calling Wildtype
}
\description{
  Decides whether a position is variant, wildtype, or uncallable,
  according to the estimated power of the given calling filters.
}
\usage{
callWildtype(reads, variants, calling.filters, pos = NULL, ...)
minCallableCoverage(calling.filters, power = 0.80, max.coverage = 1000L)
}
\arguments{
  \item{reads}{
    The read alignments, i.e., a path to a BAM file, or the coverage,
    including a \code{BigWigFile} object.
  }
  \item{variants}{
    The called variants, a tally \code{GRanges}.
  }
  \item{calling.filters}{
    Filters used to call the variants.
  }
  \item{pos}{
    A \code{GRanges} indicating positions to query; output is in the
    same order. If this is \code{NULL}, the entire genome is
    considered. This is not called \code{which}, because we are
    indicating positions, not selecting from regions.
  }
  \item{power}{
    The chance of detecting a variant if one is there.
  }
  \item{max.coverage}{
    The max coverage to be considered for the minimum (should not need
    to be tweaked).
  }
  \item{\dots}{
    Arguments to pass down to \code{minCallableCoverage}.
  }
}
\details{
  For each position (in the genome, or as specified by \code{pos}), the
  coverage is compared against the return value of
  \code{minCallableCoverage}. If the coverage is above the callable
  minimum, the position is called, either as a variant (if it is in
  \code{variants}) or wildtype. Otherwise, it is considered a no-call.

  The \code{minCallableCoverage} function
  expects and only considers the filters returned by
  \code{\link{VariantCallingFilters}}.
}
\value{
  A logical vector (or logical \code{RleList} if \code{pos} is
  \code{NULL}), that is \code{TRUE} for wildtype, \code{FALSE} for variant,
  \code{NA} for no-call.
}
\author{
  Michael Lawrence
}
\examples{
p53 <- gmapR:::exonsOnTP53Genome("TP53")
bams <- LungCancerLines::LungCancerBamFiles()
bam <- bams$H1993
tally.param <- TallyVariantsParam(gmapR::TP53Genome(), 
                                  high_base_quality = 23L,
                                  which = range(p53))
called.variants <- callVariants(bam, tally.param)

pos <- c(called.variants, shift(called.variants, 3))
wildtype <- callWildtype(bam, called.variants, VariantCallingFilters(), 
                         pos = pos, power = 0.85)
}

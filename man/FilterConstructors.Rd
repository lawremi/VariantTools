\name{FilterConstructors}
\alias{SetdiffVariantsFilter}
\alias{MinTotalDepthFilter}
\alias{MaxControlFreqFilter}
\alias{DepthFETFilter}

\title{Variant Filter Constructors}

\description{
  These functions construct filters (implemented as functions) suitable
  for collection into \code{FilterRules} objects, which are then used to
  filter variant calls. See examples.
}

\usage{
SetdiffVariantsFilter(other)
MinTotalDepthFilter(min.depth = 10L)
MaxControlFreqFilter(control, control.cov, max.control.freq = 0.03)
DepthFETFilter(control, control.cov, p.value.cutoff = 0.05)
}

\arguments{
  \item{other}{The set of variants (as a \code{VRanges}) to
    subtract from the set being filtered.
  }
  \item{min.depth}{The minimum depth for a variant to pass.}
  \item{control}{The control set of variants (as a \code{VRanges}) to
    use when filtering for case-specific variants.
  }
  \item{control.cov}{The coverage (as an \code{RleList}) for the sample
    corresponding to the calls in \code{control}.
  }
  \item{max.control.freq}{The maximum alt frequency allowed in the
    control for a variant to be considered case-specific.
  }
  \item{p.value.cutoff}{Passing variants must have a p-value below this
    value.}
}

\value{
  In all cases, a closure that returns a logical vector indicating which
  elements of its argument should be retained.
}

\author{
  Michael Lawrence
}

\seealso{
  There are some convenience functions that construct \code{FilterRules}
  objects that contain one or more of these filters. Examples are
  \code{\link{VariantQAFilters}} and
  \code{\link{VariantCallingFilters}}.
}

\examples{
  ## Find case-specific variants in a case/control study
  
  bams <- LungCancerLines::LungCancerBamFiles()
  tally.param <- TallyVariantsParam(gmapR::TP53Genome(), 
                                      high_base_quality = 23L,
                                      which = gmapR::TP53Which())

  case <- callVariants(bams$H1993, tally.param)
  control <- callVariants(bams$H2073, tally.param)
  control.cov <- coverage(bams$H2073)

  filters <-
    FilterRules(list(caseOnly = SetdiffVariantsFilter(control),
                       minTotalDepth = MinTotalDepthFilter(min.depth=10L),
                       maxControlFreq = MaxControlFreqFilter(control,
                           control.cov, max.control.freq=0.03),
                       depthFET = DepthFETFilter(control, control.cov,
                           p.value.cutoff=0.05)
                       ))
  
  specific <- subsetByFilter(case, filters)
}

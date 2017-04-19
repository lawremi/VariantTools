\name{annotateWithControlDepth}
\alias{annotateWithControlDepth}
\title{
  Annotate Case with Control Depth
}
\description{
  Matches the \code{case} variants to the (typically unfiltered)
  \code{control} variants and returns \code{case}, with the additional
  metadata columns \code{control.alt.depth} and
  \code{control.total.depth}, corresponding to \code{altDepth(control)}
  and \code{totalDepth(control)}, respectively.
}
\usage{
annotateWithControlDepth(case, control, control.cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{case}{
    The variants of interest, as a \code{VRanges}.
  }
  \item{control}{
    The control variants, typically unfiltered, as a \code{VRanges}.
  }
  \item{control.cov}{
    The control coverage, as an \code{RleList}.
  }
}
\value{
  \code{case}, plus two new metadata columns, \code{control.alt.depth}
  and \code{control.total.depth}
}
\details{
  If a \code{case} variant is not found in \code{control}, a count of
  \code{0} is assigned to \code{control.alt.depth}, under the assumption
  that the \code{control} object is not filtered, i.e., it contains the
  raw tallies.
}
\author{
  Michael Lawrence
}
\seealso{
  \code{\link{callSampleSpecificVariants}}, which uses this function.
}
\examples{
  bams <- LungCancerLines::LungCancerBamFiles()
  data(vignette)
  case <- callVariants(tallies_H1993)
  control <- tallies_H2073
  control.cov <- coverage(bams$H2073)
  annotateWithControlDepth(case, control, control.cov)
}

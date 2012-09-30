\name{callableFraction}
\alias{callCallable}
\alias{callableFraction}
\title{
  Callable fractions of transcripts 
}
\description{
  Given a set of alignments and a set of genomic features,
  \code{callableFraction} will return, for each feature, the callable
  fraction of its width. This is according to the likelihood ratio test,
  with power \code{power}. The helper function \code{callCallable}
  returns a logical \code{RleList} indicating the callable positions.
}
\usage{
callCallable(reads, plower = 0.2, perror = 1/1000, power = 0.999)
callableFraction(tx, reads, ...)
}
\arguments{
  \item{tx}{
    A \code{GRangesList} object, like one returned
    by \code{\link[GenomicFeatures]{exonsBy}} or
    \code{\link[GenomicFeatures]{cdsBy}}.
  }
  \item{reads}{
    A BAM file or \code{GappedAlignment(Pair)s} object
  }
  \item{plower}{
    The lower bound on the alternative frequency of interest.
    Should have the same value as used in the \code{\link{callVariants}}
    function.
  }
  \item{perror}{
    Assumed error rate from sequencing, should match that used by
    \code{\link{callVariants}}.
  }
  \item{prob}{
    The prefered power level for the homozygous reference
    detection. This parameter will affect the global false positive rate
    of hom-ref calls.
  }
  \item{\dots}{Arguments passed down to \code{callCallable}}
}
\value{
  Returns a named numeric vector the length of the transcripts in the txdb
  with the fraction of the transcript that has coverage above the
  callable threshold.
}

\author{
  Michael Lawrence, Jeremiah Degenhardt
}

\keyword{internal}

\name{callWildtype}
\alias{callWildtype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Call wildtype}
\description{
  A function to generate a logical RLE of homozygous ref (TRUE), variant
  (FALSE) and no call (NA) for a sample across the complete genome.
}
\usage{
callWildtype(cov, var_gr, power = 0.999, plower = 0.2, perror = 1/1000)
}
\arguments{
  \item{cov}{
    RLE vector of the coverage 
  }
  \item{var_gr}{
    Called variants to be used in generating the logical RLE
  }
  \item{power}{
    The minimum power need to be able to call a wildtype position
  }
  \item{plower}{
    The prefered alternate frequency to test against
  }
  \item{perror}{
    The assumed sequencing error rate
  }
}
\value{
  A logical RLE of called positions, TRUE = wildtype, FALSE = variant, and
  no-call positions (NA).
}
\author{
  Michael Lawrence, Jeremiah Degenhardt
}

\name{extractCoverageForPositions}
\alias{extractCoverageForPositions}
\title{
  Get Coverage at Positions
}
\description{
  Gets values from an \code{RleList} corresponding to positions (width 1
  ranges) in a \code{GRanges} (or \code{VRanges}). The result is a
  simple atomic vector.
}
\usage{
extractCoverageForPositions(cov, pos)
}
\arguments{
  \item{cov}{
    An \code{RleList} like that returned by \code{\link[IRanges]{coverage}}.
  }
  \item{pos}{
    A \code{GRanges} consisting only of width-1 ranges.
  }
}
\value{
  Atomic vector with one value from \code{cov} per position in \code{pos}.
}
\author{
  Michael Lawrence
}

\name{epGRanges}
\alias{epGRanges}
\title{
  Equipartition a GRanges object. 
}
\description{
Splits a GRanges object into a specified number of partitions with
equal (total) width. This is useful to ensure balanced
worker load in parallel evaluation of tallies for instance.
Each returned partition is a GRanges and the set of all partitions
is returned as a GRangesList object. 
}
\usage{
epGRanges(gr,partitionNb)
}
\arguments{
  \item{gr}{
    A \code{GRanges} object to be partitioned.
  }
  \item{partitionNb}{
    The desired number of partitions.
  }
}
\value{
    A \code{GRangesList} object where all elements (\code{GRanges})
    have the same (total) width.
}
\author{
    Alexandre Kuhn
}
\examples{
   gr<-GRanges(c("chr1","chr2"),IRanges(c(1,1),c(100,1e5)))
   epGRanges(gr,2)   
}

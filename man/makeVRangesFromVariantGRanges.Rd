\name{makeVRangesFromVariantGRanges}
\alias{makeVRangesFromVariantGRanges}
\title{
  Updating to VRanges
}
\description{
  This function converts a variant-style GRanges (as output by
  VariantTools 1.2.x) to a VRanges object.
}
\usage{
makeVRangesFromVariantGRanges(x, genome)
}
\arguments{
  \item{x}{
    A variant GRanges
  }
  \item{genome}{
    A GmapGenome object (for retrieving reference anchors for indels)
  }
}
\value{
  A VRanges
}
\author{
  Michael Lawrence
}

\name{plotTiTv}
\alias{plotTiTv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A function to make transition/transversion mosaic plots from variant calls
}
\description{
  This function will take a variant call Granges and produce a Ti/Tv
  mosaic plot
}
\usage{
plotTiTv(gr, main = "TiTv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr}{
    Variant call Granges as produced by callVariantsP, or tumorNormalCompare
  }
  \item{main}{
    main title passed down to plot
    }
}
\value{
  Generates a plot of the data
}
\author{
  Jeremiah Degenhardt
}

\keyword{internal}

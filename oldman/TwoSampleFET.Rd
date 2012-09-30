\name{TwoSampleFET}
\alias{TwoSampleFET}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function to find sample-specific mutations by a two-sample comparison
using the Fisher's exact test and a simple frequency threshold
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  This function takes and GRanges of called varaints from sample1 , a GRanges of raw
  variants from sample2, the coverage object from sample2 and two
  parameters and returns a GRanges of the variants specific to sample1 
}
\usage{
TwoSampleFET(gr, raw, cov, pval = 0.001, freq = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr}{
    %%     ~~Describe \code{gr} here~~
    GRanges of called variants from sample1
}
  \item{raw}{
    %%     ~~Describe \code{raw} here~~
    GRnages of raw variants from sample 2
}
  \item{cov}{
    %%     ~~Describe \code{cov} here~~
    RLE vector of coverage from sample2
}
  \item{pval}{
    %%     ~~Describe \code{pval} here~~
    the FET p-value to be used for filtering
}
  \item{freq}{
    %%     ~~Describe \code{freq} here~~
    The minimum frequency threshold needed in sample2 before the variant
    is considered to be seen
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Variant GRanges of the sample-specific variants
}
\author{
%%  ~~who you are~~
Jeremiah Degenhardt
}
\keyword{internal}

\name{checkVariantConcordance}
\alias{checkVariantConcordance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to check the variant call concordance between two samples or
lane. 
}
\description{
This function takes two variant call GRanges objects (as produced by
callVariants) and will check the concordance of the variant calls. Uses
a subset of the variants with the best coverage.
}

\usage{
checkVariantConcordance(gr1, gr2, SNP_gr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr1}{
    %%     ~~Describe \code{gr1} here~~
    Variant GRanges object from sample/lane 1 (as produced by callVariants)
  }
  \item{gr2}{
    %%     ~~Describe \code{gr2} here~~
    Variant GRanges object from sample/lane 2 (as produced by callVariants)
  }
  \item{SNP_gr}{
    a GRanges object of the SNPs you would like to constrain the check
  by
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A numeric vector with the concordance and the number of positions checked
}
\author{
  %%  ~~who you are~~
  Jeremiah Degenhardt
}
\keyword{internal}
% \examples{
% }

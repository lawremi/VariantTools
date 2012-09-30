\name{getChrGRangesFromGmap}
\alias{getChrGRangesFromGmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function to retreive a GRanges object of the chromosome names, start
and end positions for the genome being used by gSNAP
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given the gSNAP genome and genome directory information return a GRanges
of the chromosome information
}
\usage{
getChrGRangesFromGmap(genome, genome_dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genome}{
    %%     ~~Describe \code{genome} here~~
    the name of the gSNAP genome that you want to get chromosome
    information from
  }
  \item{genome_dir}{
    %%     ~~Describe \code{genome_dir} here~~
    The directory containing the gSNAP genome 
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
GRanges object with the chromosome names, start and end position
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  %%  ~~who you are~~
  Jeremiah Degenhardt
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{vignette}
\alias{tallies_H1993}
\alias{tallies_H2073}
\alias{coverage_H1993}
\alias{coverage_H2073}
\alias{genome_p53}
\alias{p53}

\docType{data}
\title{Vignette Data}
\description{
  Precomputed data for use in the vignette, mostly for the sake of
  Windows, where gmapR and its tallying functionality are unsupported.
}
\usage{data(vignette)}
\format{
  The following objects are included:
  \describe{
    \item{tallies_H1993, tallies_H2073}{Tallies for the two samples.}
    \item{coverage_H1993, coverage_H2073}{Coverage for the two samples.}
    \item{p53}{A GRanges of the p53 exons}
    \item{genome_p53}{DNAStringSet with the genome sequence of the p53
      region}
  }
}
\details{
  The following demonstrates how we created these objects:
  \preformatted{
    bams <- LungCancerLines::LungCancerBamFiles()
    tally.param <- TallyVariantsParam(gmapR::TP53Genome(), 
                                      high_base_quality = 23L,
				      which = range(p53) + 5e4,
                                      indels = TRUE, read_length = 75L)
    tallies_H1993 <- tallyVariants(bams$H1993, tally.param)
    tallies_H2073 <- tallyVariants(bams$H2073, tally.param)
    coverage_H1993 <- coverage(bams$H1993)
    coverage_H2073 <- coverage(bams$H2073)
    genome_p53 <- DNAStringSet(getSeq(gmapR::TP53Genome()))
    p53 <- gmapR:::exonsOnTP53Genome("TP53")
  }
}
\source{
  Computed from the data in the LungCancerLines package.
}
\examples{
data(vignette)
}
\keyword{datasets}

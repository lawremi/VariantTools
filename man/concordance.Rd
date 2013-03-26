\name{concordance}
\alias{calculateConcordanceMatrix}
\alias{calculateVariantConcordance}
\alias{callVariantConcordance}
\title{
  Variant Concordance
}
\description{
  Functions for calculating concordance between variant sets and
  deciding whether two samples have identical genomes.
}
\usage{
calculateVariantConcordance(gr1, gr2, which = NULL)
calculateConcordanceMatrix(variantFiles)
callVariantConcordance(concordanceMatrix, threshold)
}
\arguments{
  \item{gr1, gr2}{The two tally \code{GRanges} to compare}
  \item{which}{
    A \code{GRanges} of positions to which the comparison is limited.
  }
  \item{variantFiles}{
    Character vector of paths to serialized tally \code{GRanges}.
  }
  \item{concordanceMatrix}{
    A matrix of concordance fractions between sample pairs, as returend by
    \code{calculateConcordanceMatrix}.
  }
  \item{threshold}{
    The concordance fraction above which edges are generated between
    samples when forming the graph.
  }
}
\details{
  The \code{calculateVariantConcordance} calculates the fraction of
  concordant variants between two samples. Concordance is defined as
  having the same position and alt allele.

  The \code{calculateConcordanceMatrix} function generates a numeric
  matrix with the concordance for each pair of samples. It accepts paths
  to serialized objects so that all variant calls are not loaded in
  memory at once. This probably should support VCF files, eventually.

  The \code{callVariantConcordance} function generates a
  concordant/non-concordant/undecidable status for each sample (that are
  assumed to originate from the same individual), given the output of
  \code{calculateConcordanceMatrix}. The status is decided as follows. A
  graph is formed from the concordance matrix using \code{threshold} to
  generate the edges. If there are multiple cliques in the graph that
  each have more than one sample, every sample is declared
  undecidable. Otherwise, the samples in the clique with more than one
  sample, if any, are marked as concordant, and the others (in singleton
  cliques) are marked as discordant.
}
\value{
  Fraction of concordant variants for \code{calculateVariantConcordance}, a
  numeric matrix of concordances for \code{calculateConcordanceMatrix},
  or a character vector of status codes, named by sample, for
  \code{callVariantConcordance}.
}
\author{
  Cory Barr (code), Michael Lawrence (inferred documentation)
}

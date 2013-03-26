\name{matchVariants}
\alias{matchVariants}
\alias{\%variant_in\%}
\title{
  Match variants by position and allele
}
\description{
  This function behaves like \code{match}, where two elements match when
  they share the same position and \dQuote{alt} allele.
}
\usage{
matchVariants(x, table)
x \%variant_in\% table
}
\arguments{
  \item{x}{
    The variants (\code{GRanges}) to match into \code{table}; the alt
    allele must be in the \dQuote{alt} metacolumn.
  }
  \item{table}{
    The variants (\code{GRanges}) to be matched into; the alt
    allele must be in the \dQuote{alt} metacolumn.
  }
}
\value{
  For \code{matchVariants}, an integer vector with the matching index in
  \code{table} for each variant in \code{x}, or \code{NA} if there is no
  match. For \code{\%variant_in\%}, a logical vector indicating whether
  there was such a match.
}
\author{
  Michael Lawrence
}

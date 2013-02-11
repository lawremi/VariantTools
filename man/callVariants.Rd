\name{callVariants}
\alias{callVariants}
\alias{callVariants,BamFile-method}
\alias{callVariants,character-method}
\alias{callVariants,GenomicRanges-method}
\alias{VariantCallingFilters}

\title{
  Call Variants
}
\description{
  Calls variants from either a BAM file or a tally \code{GRanges}
  object. The variants are first filtered with \code{\link{qaVariants}},
  and the remaining candidates are called using a binomial likelihood
  ratio test: P(D|p=\code{p.lower}) / P(D|p=\code{p.error}) > 1.
}
\usage{
\S4method{callVariants}{BamFile}(x, tally.param,
                                 qa.filters = VariantQAFilters(),
                                 calling.filters = VariantCallingFilters(...),
                                 ...)
\S4method{callVariants}{character}(x, ...)
\S4method{callVariants}{GenomicRanges}(x,
     calling.filters = VariantCallingFilters(...), 
     ...)
VariantCallingFilters(read.count = 2L, p.lower = 0.2, p.error = 1e-4)
}

\arguments{
  \item{x}{
    Either a path to an indexed bam, a \code{BamFile}
    object, or a \code{GRanges} as returned by \code{\link{tallyVariants}}. 
  }
  \item{tally.param}{
    Parameters controlling the variant tallying step,
    as typically constructed by \code{\link{VariantTallyParam}}.
  }
  \item{qa.filters}{
    Filters used in the QA step, see \code{\link{VariantQAFilters}}.
  }
  \item{calling.filters}{
    Filters used in the calling step, typically constructed with
    \code{VariantCallingFilters}, see arguments listed below.
  }
  \code{...}{Arguments for \code{VariantCallingFilters}, listed below.}
  \item{read.count}{
    Require at least this many high quality reads with
    the alternate base. The default value is designed to catch sequencing
    errors where coverage is too low to rely on the LRT. Increasing this
    value has a significant negative impact on power.
  }
  \item{p.lower}{The lower bound on the binomial probability for a true
    variant.
  }
  \item{p.error}{The binomial probability for a sequencing error (default
    is reasonable for Illumina data with the default quality cutoff).
  }
  \item{\dots}{
    Arguments to pass to \code{VariantCallingFilters}.
  }
}
\value{
  For \code{callVariants}, a \code{GRanges} of the called variants (the
  tallies that pass the QA and calling filters). See the documentation
  of \code{\link[gmapR]{bam_tally}} for complete details.
  
  For \code{VariantCallingFilters}, a \code{\link[IRanges]{FilterRules}}
  object with the filters for calling the variants (presumably after the
  QA filters have been applied).
}
\author{
  Michael Lawrence, Jeremiah Degenhardt
}
\examples{
bams <- LungCancerLines::LungCancerBamFiles()
tally.param <- VariantTallyParam(gmapR::TP53Genome(), 
                                 readlen = 100L,
                                 high_base_quality = 23L,
                                 which = gmapR::TP53Which())

## simple usage
variants <- callVariants(bams$H1993, tally.param)

## customize
qa.filters <- VariantQAFilters(fisher.strand.p.value = 1e-4)
calling.filters <- VariantCallingFilters(p.error = 1/1000)
callVariants(bams$H1993, tally.param, qa.filters, calling.filters)
}

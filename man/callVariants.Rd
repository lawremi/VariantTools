\name{callVariants}
\alias{callVariants}
\alias{callVariants,BamFile-method}
\alias{callVariants,character-method}
\alias{callVariants,GenomicRanges-method}
\alias{callVariants,VRanges-method}
\alias{VariantCallingFilters}

\title{
  Call Variants
}
\description{
  Calls variants from either a BAM file or a \code{VRanges}
  object. The variants are called using a binomial likelihood
  ratio test. Those calls are then subjected to a post-filtering step.
}
\usage{
\S4method{callVariants}{BamFile}(x, tally.param,
                                 calling.filters = VariantCallingFilters(...),
                                 post.filters = FilterRules(),
                                 ...)
\S4method{callVariants}{character}(x, ...)
\S4method{callVariants}{VRanges}(x,
     calling.filters = VariantCallingFilters(...),
     post.filters = FilterRules(),
     ...)
VariantCallingFilters(read.count = 2L, p.lower = 0.2, p.error = 1/1000)
}

\arguments{
  \item{x}{
    Either a path to an indexed bam, a \code{BamFile}
    object, or a \code{VRanges} as returned by \code{\link{tallyVariants}}. 
  }
  \item{tally.param}{
    Parameters controlling the variant tallying step,
    as typically constructed by \code{\link{TallyVariantsParam}}.
  }
  \item{calling.filters}{
    Filters used in the calling step, typically constructed with
    \code{VariantCallingFilters}, see arguments listed below.
  }
  \item{post.filters}{
    Filters that are applied after the initial calling step. These
    consider the set of variant calls as a whole and remove those with
    suspicious patterns.
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
\details{
  There are two steps for calling
  variants: the actual statistical test that decides whether a variant
  exists in the data, and a post-filtering step. By default, the initial
  calling is based on a binomial likelihood ratio test
  (P(D|p=\code{p.lower}) / P(D|p=\code{p.error}) > 1). The test amounts
  to excluding putative variants with less than ~4\% alt frequency. A
  variant is also required to be represented by at least 2 alt
  reads. The post-filtering stage considers the set of variant calls as
  a whole and removes variants with suspicious patterns. Currently,
  there is a single post-filter, disabled by default, that removes
  variants that are clumped together on the chromosome (see the
  \code{max.nbor.count} parameter).
}
\value{
  For \code{callVariants}, a \code{VRanges} of the called variants (the
  tallies that pass the calling filters). See the documentation
  of \code{\link[gmapR]{bam_tally}} for complete details.
  
  For \code{VariantCallingFilters}, a \code{\link[IRanges]{FilterRules}}
  object with the filters for calling the variants.
}
\author{
  Michael Lawrence, Jeremiah Degenhardt
}
\examples{
bams <- LungCancerLines::LungCancerBamFiles()
tally.param <- TallyVariantsParam(gmapR::TP53Genome(), 
                                  high_base_quality = 23L,
                                  which = gmapR::TP53Which())

## simple usage
variants <- callVariants(bams$H1993, tally.param)

## customize
calling.filters <- VariantCallingFilters(p.error = 1/1000)
callVariants(bams$H1993, tally.param, calling.filters)
}

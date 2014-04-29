\name{callSampleSpecificVariants}
\alias{callSampleSpecificVariants}
\alias{callSampleSpecificVariants,BamFile,BamFile-method}
\alias{callSampleSpecificVariants,character,character-method}
\alias{callSampleSpecificVariants,GenomicRanges,GenomicRanges-method}
\alias{callSampleSpecificVariants,VRanges,VRanges-method}
\alias{SampleSpecificVariantFilters}

\title{Call Sample-Specific Variants}
\description{
  Calls sample-specific variants by comparing case and control variants
  from paired samples, starting from the BAM files or unfiltered
  tallies. For example, these variants would be considered somatic
  mutations in a tumor vs. normal comparison.
}
\usage{
SampleSpecificVariantFilters(control, control.cov, calling.filters,
                             power = 0.8, p.value = 0.01)
\S4method{callSampleSpecificVariants}{BamFile,BamFile}(case, control,
  tally.param, calling.filters = VariantCallingFilters(), post.filters =
  FilterRules(), ...)
\S4method{callSampleSpecificVariants}{character,character}(case, control, ...)
\S4method{callSampleSpecificVariants}{VRanges,VRanges}(case,
  control, control.cov, ...)
## DEPRECATED
\S4method{callSampleSpecificVariants}{GenomicRanges,GenomicRanges}(case,
  control, control.cov,
  calling.filters = VariantCallingFilters(), post.filters =
  FilterRules(), ...)
}
\arguments{
  \item{case}{
    The BAM file for the case, or the called variants as output by
    \code{\link{callVariants}}.
  }
  \item{control}{
    The BAM file for the control, or the raw tallies as output by
    \code{\link{tallyVariants}}.
  }
  \item{tally.param}{Parameters controlling the variant tallying step,
    as typically constructed by \code{\link{TallyVariantsParam}}.}
  \item{calling.filters}{Filters to use for the initial,
    single-sample calling against reference, typically constructed by
    \code{\link{VariantCallingFilters}}.
  }
  \item{post.filters}{
    Filters that are applied after the initial calling step. These
    consider the set of variant calls as a whole and remove those with
    suspicious patterns. They are only applied to the \code{case}
    sample; only QA filters are applied to \code{control}.
  }
  \item{\dots}{
    For a BAM file, arguments to pass down to the \code{GenomicRanges}
    method. For the \code{GenomicRanges} method, arguments to pass down
    to \code{SampleSpecificVariantFilters}, except for
    \code{control.cov}, \code{control.called}, \code{control.raw} and
    \code{lr.filter}.
  }
  \item{control.cov}{The coverage for the control sample.}
  \item{power}{The power cutoff, beneath which a variant will not be
    called case-specific, due to lack of power in control.
  }
  \item{p.value}{The binomial p-value cutoff for determining whether the
    control frequency is sufficiently extreme (low) compared to the case
    frequency. A p-value below this cutoff means that the variant will
    be called case-specific.
  }
}
\details{
  For each sample, the variants are tallied (when the input is BAM), QA
  filtered (case only), called and determined to be sample-specific.
  The \code{callSampleSpecificVariants} function is fairly high-level,
  but it still allows the user to override the parameters and filters
  for each stage of the process. See \code{\link{TallyVariantsParam}},
  \code{\link{VariantQAFilters}}, \code{\link{VariantCallingFilters}}
  and \code{SampleSpecificVariantFilters}.

  It is safest to pass a BAM file, so that the computations are
  consistent for both samples. The \code{GenomicRanges} method is
  provided mostly for optimization purposes, since tallying the variants
  over the entire genome is time-consuming. For small gene-size regions,
  performance should not be a concern.
  
  This is the algorithm that determines whether a variant is specific to
  the case sample:
  \enumerate{
    \item{Filter out all case calls that were also called in
      control. The \code{callSampleSpecificVariants} function does
      \strong{not} apply the QA filters when calling variants in
      control. This prevents a variant from being called specific to
      case merely due to questionable data in the control.
    }
    \item{For the remaining case calls, calculate whether there was
      sufficient power in control under the likelihood ratio test, for a
      variant present at the \code{p.lower} frequency. If that is below
      the \code{power} cutoff, discard it.
    }
    \item{For the remaining case calls, test whether the control
      frequency is sufficient extreme (low) compared to the case
      frequency, under the binomial model. The null hypothesis is that
      the frequencies are the same, so if the test p-value is above
      \code{p.value}, discard the variant. Otherwise, the variant is
      called case-specific.
    }
  }
}
\value{
  A \code{VRanges} with the case-specific variants (such as
  somatic mutations).
}
\author{
  Michael Lawrence, Jeremiah Degenhardt
}
\examples{
bams <- LungCancerLines::LungCancerBamFiles()
tally.param <- TallyVariantsParam(gmapR::TP53Genome(), 
                                  high_base_quality = 23L,
                                  which = gmapR::TP53Which())
callSampleSpecificVariants(bams$H1993, bams$H2073, tally.param)
}

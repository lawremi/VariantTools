\name{callGenotypes}
\alias{callGenotypes}
\alias{callGenotypes,VRanges,BigWigFile-method}
\alias{callGenotypes,TabixFile,BigWigFile-method}
\alias{CallGenotypesParam}

\title{
  Call Genotypes
}
\description{
  Calls genotypes from a set of tallies (such as a \code{VRanges} or VCF
  file) and the coverage (currently as a \code{BigWigFile}). We call the
  genotype with the highest likelihood, where the likelihood is based on
  a binomial model of the variant frequency.
}
\usage{
\S4method{callGenotypes}{VRanges,BigWigFile}(variants, cov,
  param = CallGenotypesParam(GmapGenome(unique(genome(variants)))),
  BPPARAM = defaultBPPARAM())
\S4method{callGenotypes}{TabixFile,BigWigFile}(variants, cov,
  param = CallGenotypesParam(GmapGenome(unique(genome(variants)))),
  BPPARAM = defaultBPPARAM())
CallGenotypesParam(genome,
                   gq.breaks = c(0, 5, 20, 60, Inf),
                   p.error = 0.05,
                   which = tileGenome(seqinfo(genome), ntile=ntile),
                   ntile = 100L)
}

\arguments{
  \item{variants}{
    Either \code{VRanges} as returned by \code{\link{tallyVariants}}, or
    a \code{TabixFile} object pointing to a VCF file. Typically, these
    tallies are \emph{not} filtered by e.g. \code{\link{callVariants}},
    because it would seem more appropriate to filter on the genotype
    quality.
  }
  \item{cov}{
    The coverage, stored as a \code{BigWigFile}.
  }
  \item{param}{
    Parameters controlling the genotyping, constructed by
    \code{\link{CallGenotypesParam}}. The default value uses the genome
    from \code{variants}.
  }
  \item{genome}{
    The \code{GmapGenome} object used during tallying.
  }
  \item{gq.breaks}{A numeric vector representing an increasing sequence
    of genotype quality breaks to segment the wildtype runs.
  }
  \item{p.error}{The binomial probability for an error. This is used to
    calculate the expected frequency of hom-ref and hom-alt variants.
  }
  \item{which}{A \code{GenomicRangesList} indicating the genomic regions
    in which to compute the genotypes. The default is to partition the
    \code{genome} into \code{ntile} tiles.
  }
  \item{ntile}{When \code{which} is missing, this indicates the number
    of tiles to generate from the genome.
  }
  \item{BPPARAM}{A \code{BiocParallelParam} object communicating the
    parallelization strategy. One job is created per tile.}
}
\details{
  In general, the behavior is very similar to that of the GATK
  UnifiedGenotyper (see references). For every position in the tallies,
  we compute a binomial likelihood for each of wildtype (0/0), het (0/1)
  and hom-alt (1/1), assuming the alt allele frequency to be
  \code{p.error}, \code{0.5} and \code{1 - p.error}, respectively. The
  genotype with the maximum likelihood is chosen as the genotype, and
  the genotype quality is computed by taking the fraction of the maximum
  likelihood over the sum of the three likelihoods.

  We assume that any position not present in the input tallies is
  wildtype (0/0) and compute the quality for every such position, using
  the provided coverage data. For scalability reasons, we segment runs
  of these positions according to user-specified breaks on the genotype
  quality. The segments become special records in the returned
  \code{VRanges}, where the range represents the segment, the \code{ref}
  is the first reference base, \code{alt} is \code{<NON_REF>} and the
  \code{totalDepth} is the mean of the coverage.

  The genotype information is recorded as metadata columns named
  according to gVCF conventions:
  \describe{
    \item{\code{GT}}{The genotype call string: \code{0/0}, \code{0/1},
      \code{1/1}.
    }
    \item{\code{GQ}}{The numeric genotype quality, phred scaled. For
      wildtype runs, this is minimum over the run.
    }
    \item{\code{PL}}{A 3 column matrix with the likelihood for each
      genotype, phred scaled. We take the minimum over wildtype runs.
    }
    \item{\code{MIN_DP}}{The minimum coverage over a wildtype
      run; \code{NA} for single positions.
    }
  }  
}
\value{
  For \code{callGenotypes}, a \code{VRanges} annotated with the genotype
  call information, as described in the details.
}
\author{
  Michael Lawrence
}
\references{
  The Genome Analysis Toolkit: a MapReduce framework for analyzing
  next-generation DNA sequencing data McKenna A, Hanna M, Banks E,
  Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D,
  Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303.
}
\examples{
## generate tallies
library(VariantTools)
bams <- LungCancerLines::LungCancerBamFiles()
tally.param <- TallyVariantsParam(gmapR::TP53Genome(), 
                                  high_base_quality = 23L,
                                  which = gmapR::TP53Which())
tallies <- tallyVariants(bams$H1993, tally.param)
sampleNames(tallies) <- "H1993"
mcols(tallies) <- NULL

## we assume the coverage is already stored as a BigWigFile
bwf <- rtracklayer::export(coverage(bams$H1993), "coverage.bw")

## simple usage
genotypes <- callGenotypes(tallies, bwf, BPPARAM=BiocParallel::SerialParam())

## write to gVCF
writeVcf(genotypes, "genotypes.vcf.gz")

## customize
params <- CallGenotypesParam(gmapR::TP53Genome(), p.error = 1/1000)
callGenotypes(tallies, bwf, params)
}

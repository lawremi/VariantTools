\name{tallyVariants}
\alias{tallyVariants}
\alias{tallyVariants,character-method}
\alias{tallyVariants,BamFile-method}
\alias{VariantTallyParam}

\title{
  Tally the positions in a BAM file
}
\description{
  Tallies the bases, qualities and read positions for every genomic
  position in a BAM file. By default, this only returns the positions
  for which an alternate base has been detected.  The typical usage is
  to pass a BAM file, the genome, the (fixed) \code{readlen} and (if the
  variant calling should consider quality) an appropriate
  \code{high_base_quality} cutoff. Passing a \code{which} argument
  allows computing on only a subregion of the genome.
}
\usage{
\S4method{tallyVariants}{BamFile}(x, param = VariantTallyParam(...), ...,
                                  BPPARAM = defaultBPPARAM())
\S4method{tallyVariants}{character}(x, ...)
VariantTallyParam(genome, readlen = NA,
                  read_pos_flank_width = 10L,
                  read_pos_breaks = flankingCycleBreaks(readlen,
                  read_pos_flank_width),
                  high_base_quality = 0L,
                  minimum_mapq = 13L,
                  variant_strand = 1L, ignore_query_Ns = TRUE,
                  ignore_duplicates = TRUE,
                  ...)
}
\arguments{
  \item{x}{
    An indexed BAM file, either a path or a \code{BamFile} object.
  }
  \item{param}{
    The parameters for the tallying process, as a
    \code{\link[gmapR]{BamTallyParam}}, typically constructed with
    \code{VariantTallyParam}, see arguments below.
  }
  \item{\dots}{
    For \code{tallyVariants}, arguments to pass to
    \code{VariantTallyParam}, listed below. For
    \code{VariantTallyParam}, arguments to pass to
    \code{\link[gmapR]{BamTallyParam}}.
  }
  \item{genome}{The genome, either a \code{\link[gmapR]{GmapGenome}} or
    something coercible to one.
  }
  \item{readlen, read_pos_flank_width}{
    If \code{read_pos_breaks} is missing, these two arguments are used to
    generate a \code{read_pos_breaks} for three bins, with the two outside
    bins having \code{read_pos_flank_width}. If \code{readlen} is
    \code{NA}, \code{read_pos_breaks} is not generated.
  }
  \item{read_pos_breaks}{The breaks used for tabulating the read positions (read
    positions) at each position. If this information is included (not
    \code{NULL}), \code{\link{qaVariants}} will use it during filtering.
  }
  \item{high_base_quality}{The minimum cutoff for whether a base is
    counted as high quality. By default, \code{\link{callVariants}} will
    use the high quality counts in the likelihood ratio test. Note that
    \code{bam_tally} will shift your quality scores by 33 no matter what
    type they are. If Illumina (pre 1.8) this will result in a range of
    31-71. If Sanger/Illumina1.8 this will result in a range of
    0-40/41. The default counts all bases as high quality. We typically
    use 56 for old Illumina, 23 for Sanger/Illumina1.8.
  }
  \item{minimum_mapq}{Minimum MAPQ of a read for it to be included in
    the tallies. This depend on the aligner; the default is reasonable
    for \code{\link[gmapR]{gsnap}}.
  }
  \item{variant_strand}{On how many strands must an alternate base
    be detected for a position to be returned. Highly recommended to set
    this to at least 1 (otherwise, the result is huge and includes
    many uninteresting reference rows).
  }
  \item{ignore_query_Ns}{Whether to ignore N calls in the
    reads. Usually, there is no reason to set this to \code{FALSE}. If
    it is \code{FALSE}, beware of low quality datasets returning
    enormous results.
  }
  \item{ignore_duplicates}{whether to ignore reads flagged as
    PCR/optical duplicates}
  \item{BPPARAM}{A
    \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
    object specifying the resources and strategy for parallelizing the
    tally operation over the chromosomes.
  }
}
\value{
  For \code{tallyVariants}, the tally \code{GRanges}.

  For \code{VariantTallyParam}, an object with parameters suitable for
  variant calling.
}
\author{
  Michael Lawrence, Jeremiah Degenhardt
}

\examples{
tally.param <- TallyVariantsParam(gmapR::TP53Genome(), 
                                  high_base_quality = 23L,
                                  which = gmapR::TP53Which())
bams <- LungCancerLines::LungCancerBamFiles()
raw.variants <- tallyVariants(bams$H1993, tally.param)
}

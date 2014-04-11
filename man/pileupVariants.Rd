\name{pileupVariants}
\alias{pileupVariants}
\title{
  Nucleotide pileup from alignments
}
\description{
  This is an alternative to \code{\link{tallyVariants}} for generating a
  \code{VRanges} from a set of alignments (BAM file) by counting the
  nucleotides at each position. This function uses the samtools-based
  \code{\link[Rsamtools]{applyPileups}} function, instead of
  \code{\link[gmapR]{bam_tally}}. Fewer dependencies, with fewer
  statistics (none beyond the fixed columns) available in the output.
}
\usage{
pileupVariants(bams, genome, param = PileupParam(), minAltDepth = 1L,
               baseOnly = TRUE, BPPARAM = defaultBPPARAM())
}
\arguments{
  \item{bams}{
    A vector/list of BAM files as interpreted
    by \code{\link[Rsamtools]{PileupFiles}}.
  }
  \item{genome}{
    An object that provides sequence information via
    \code{\link[Biostrings]{getSeq}}.
  }
  \item{param}{
    A \code{\link[Rsamtools]{PileupParam}} object that specifies the
    mode of iteration and various filters.
  }
  \item{minAltDepth}{
    Minimal alt depth to be included in the output. The default avoids
    outputting results for positions/alleles that show no differences.
  }
  \item{baseOnly}{
    Whether to drop records with \dQuote{N} in either the ref or alt.
  }
  \item{BPPARAM}{
    Not yet supported.
  }
}
\value{
  A VRanges object with read depth information for each position,
  allele, and sample.
}
\author{
  Michael Lawrence
}
\seealso{
  \code{\link{tallyVariants}} for more statistics.
}
\examples{
     bams <- LungCancerLines::LungCancerBamFiles()
     pileup <- pileupVariants(bams, gmapR::TP53Genome(),
                              Rsamtools::PileupParam(which=gmapR::TP53Which()))
}

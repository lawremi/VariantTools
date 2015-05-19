\name{variantGR2Vcf}
\alias{variantGR2Vcf}
\title{
  Create a VCF for some variants
}
\description{
  The \strong{deprecated} way to create a
  \code{\link[VariantAnnotation]{VCF}} object from a variant/tally
  \code{GRanges}. This can then be output to a file using
  \code{\link[VariantAnnotation]{writeVcf}}. The flavor of VCF is
  specific for calling variants, not genotypes; see below.
}
\usage{
variantGR2Vcf(x, sample.id, project = NULL,
              genome = unique(GenomicRanges::genome(x)))
}
\note{
  This function is \strong{DEPRECATED}. The \code{callVariants} function
  now returns a \code{\link[VariantAnnotation]{VRanges}} object that can
  be coerced to a \code{VCF} object via \code{as(x, "VCF")}.
}
\arguments{
  \item{x}{
    The variant/tally \code{GRanges}.
  }
  \item{sample.id}{
    Unique ID for the sample in the VCF.
  }
  \item{project}{
    Description of the project/experiment; will be included in the VCF header.
  }
  \item{genome}{
    \code{GmapGenome} object, or the name of one (in the default genome
    directory). This is used for obtaining the anchor base when
    outputting indels.
  }
}
\details{
  A variant \code{GRanges} has an element for every unique combination
  of position and alternate base. A \code{VCF} object, like the file
  format, has a row for every position, with multiple alternate alleles
  collapsed within the row. This is the fundamental difference between
  the two data structures. We feel that the \code{GRanges} is easier to
  manipulate for filtering tasks, while \code{VCF} is obviously
  necessary for communication with external databases and tools.

  Normally, despite its name, VCF is used for communicating
  \emph{genotype} calls. We are calling \emph{variants}, not genotypes,
  so we have extended the format accordingly.

  Here is the mapping in detail:
  \itemize{
    \item{The \code{rowRanges} is formed by dropping the metadata columns
      from the \code{GRanges}.
    }
    \item{The \code{colData} consists of a single column,
      \dQuote{Samples}, with a single row, set to 1 and named
      \code{sample.id}.
    }
    \item{The \code{exptData} has an element \dQuote{header} with element
      \dQuote{reference} set to the \code{seqlevels(x)} and element
      \dQuote{samples} set to \code{sample.id}. This will also include the
      necessary metadata for describing our extensions to the format.
    }
    \item{The \code{fixed} table has the \dQuote{REF} and \dQuote{ALT}
      alleles, with \dQuote{QUAL} and \dQuote{FILTER} set to \code{NA}.
    }
    \item{The \code{geno} list has six matrix elements, all with a
      single column. The first is the mandatory \dQuote{GT} element, the
      genotype, which we set to \code{NA}. Then there is \dQuote{AD}
      (list matrix with the read count for each REF and ALT),
      \dQuote{DP} (integer matrix with the total read count), and
      \dQuote{AP} (list matrix of 0/1 flags for whether whether REF
      and/or ALT was present in the data).
    }
  }
}
\value{
  A \code{VCF} object.
}
\author{
  Michael Lawrence, Jeremiah Degenhardt
}
\examples{
\dontrun{
vcf <- variantGR2Vcf(variants, "H1993", "example")
writeVcf(vcf, "H1993.vcf", index = TRUE)
}
}

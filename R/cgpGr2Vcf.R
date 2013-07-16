## Convert variant GRanges into a VCF
## **DEPRECATED** in favor of VRanges

variantGR2Vcf <- function(x, sample.id, project = NULL,
                          genome = unique(GenomicRanges::genome(x)))
{
  variantGRangesIsDeprecated("variantGR2Vcf")
  vr <- variantGRangesToVRanges(x)
  sampleNames(vr) <- sample.id
  genome(vr) <- genome
  metadata(vr)$project <- project
  asVCF(vr, meta = if (!is.null(project)) "project" else character())
}


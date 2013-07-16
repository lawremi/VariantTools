### =========================================================================
### Tallying Alignments for Variant Calling
### -------------------------------------------------------------------------

setGeneric("tallyVariants", function(x, ...) standardGeneric("tallyVariants"))

setMethod("tallyVariants", "BamFile",
          function(x, param = TallyVariantsParam(...), ...,
                   mc.cores = getOption("mc.cores", 2))
          {
            tally_region <- function(which) {
              iit <- bam_tally(x, param@bamTallyParam, which = which)
              variantSummary(iit, param@read_pos_breaks,
                             param@high_base_quality)
            }
            which <- param@bamTallyParam@which
            if (length(which) == 0L) {
              genome <- param@bamTallyParam@genome
              si <- intersect(seqinfo(genome), seqinfo(x))
              ans <- applyByChromosome(si, tally_region,
                                       mc.cores = mc.cores)
            } else {
              if (length(which) == 1L) {
                chunks <- breakInChunks(width(which),
                                        ceiling(width(which) / mc.cores))
                which <- GRanges(seqnames(which),
                                 IRanges(start(which) + start(chunks) - 1L,
                                         width = width(chunks)))
              }
              ind <- seq_len(length(which))
              ans <- safe_mclapply(ind, function(i) {
                tally_region(which[i])
              }, mc.cores = mc.cores)
            }
            ans.flat <- unlist(ans, use.names = FALSE)
            ans.flat[!(ans.flat %over% param@mask)]
          })

setMethod("tallyVariants", "BamFileList", function(x, ...) {
  stackSamples(VRangesList(lapply(x, tallyVariants, ...)))
})

setMethod("tallyVariants", "character", function(x, ...) {
  tallyVariants(BamFile(x), ...)
})

## PLAN:
## - Create an IIT class that wraps the IIT externalptr
## - Rename VariantTallyParam to TallyVariantsParam (consistent with function)
## - Extract variant summary stuff from bam_tally into summarizeAlleles().
##   The bam_tally function will return an IIT object. 
## - Create an actual class named TallyVariantsParam, composed of:
##   - BamTallyParam
##   - SummarizeAllelesParam (may be overkill)
##   - Mask
## - In theory, could make tallyVariants dispatch on both the 'x' and 'param'.
##   This would allow totally new tallying algorithms.

setClass("TallyVariantsParam",
         representation(bamTallyParam = "BamTallyParam",
                        read_pos_breaks = "integer",
                        high_base_quality = "integer",
                        mask = "GenomicRanges"))

TallyVariantsParam <- function(genome,
                               read_pos_breaks = NULL,
                               high_base_quality = 0L,
                               minimum_mapq = 13L, 
                               variant_strand = 1L, ignore_query_Ns = TRUE,
                               ignore_duplicates = TRUE,
                               mask = GRanges(),
                               ...)
{
  if (variant_strand < 1 || variant_strand > 2)
    stop("'variant_strand' must be 1 or 2")
  if (!isTRUE(ignore_query_Ns))
    stop("'ignore_query_Ns' must be TRUE")
  bam.tally.args <- list(genome = genome,
                         variant_strand = variant_strand,
                         ignore_query_Ns = ignore_query_Ns,
                         minimum_mapq = minimum_mapq,
                         ignore_duplicates = ignore_duplicates,
                         ...)
  bam.tally.param <- do.call(BamTallyParam, bam.tally.args)
  new("TallyVariantsParam", bamTallyParam = bam.tally.param,
      read_pos_breaks = as.integer(read_pos_breaks),
      high_base_quality = as.integer(high_base_quality),
      mask = mask)
}

VariantTallyParam <- function(genome,
                              readlen = NA,
                              read_pos_flank_width = 10L,
                              read_pos_breaks = flankingCycleBreaks(readlen,
                                read_pos_flank_width),
                              high_base_quality = 0L,
                              minimum_mapq = 13L, 
                              variant_strand = 1L, ignore_query_Ns = TRUE,
                              ignore_duplicates = TRUE,
                              ...)
{
  .Deprecated("TallyVariantsParam")
  TallyVariantsParam(genome, read_pos_breaks,
                     high_base_quality, minimum_mapq, 
                     variant_strand, ignore_query_Ns,
                     ignore_duplicates,
                     ...)
}

setMethod("show", "TallyVariantsParam", function(object) {
  cat("A", class(object), "object\n", sep = " ")
  cat(gmapR:::showSlots(object@bamTallyParam, count = FALSE), sep = "")
  cat(gmapR:::showSlots(object, exclude = "bamTallyParam", count = FALSE),
      sep = "")
})

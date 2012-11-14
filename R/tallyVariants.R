### =========================================================================
### Tallying Alignments for Variant Calling
### -------------------------------------------------------------------------

setGeneric("tallyVariants", function(x, ...) standardGeneric("tallyVariants"))

setMethod("tallyVariants", "BamFile",
          function(x, param = VariantTallyParam(...), ...,
                   mc.cores = getOption("mc.cores", 2))
          {
            bam_tally_region <- function(which) {
              bam_tally(x, param, which = which)
            }
            which <- param@which # we depend on this being one range
            if (is.null(space(which))) {
              genome <- param@genome
              si <- intersect(seqinfo(genome), seqinfo(x))
              ans <- applyByChromosome(si, bam_tally_region,
                                       mc.cores = mc.cores)
            } else {
              which <- as(which, "GRanges")
              if (length(which) == 1L) {
                chunks <- breakInChunks(width(which),
                                        ceiling(width(which) / mc.cores))
                which <- GRanges(seqnames(which),
                                 IRanges(start(which) + start(chunks) - 1L,
                                         width = width(chunks)))
              }
              ind <- seq_len(length(which))
              ans <- safe_mclapply(ind, function(i) {
                bam_tally_region(which[i])
              }, mc.cores = mc.cores)
            }
            unlist(ans, use.names = FALSE)
          })

setMethod("tallyVariants", "character", function(x, ...) {
  tallyVariants(BamFile(x), ...)
})

## PLAN:
## - Create an IIT class that wraps the IIT externalptr
## - Rename VariantTallyParam to TallyVariantsParam (consistent with function)
## - Extract variant summary stuff from bam_tally into variantTallySummary().
##   The bam_tally function will return an IIT object. 
## - Create an actual class named TallyVariantsParam, composed of two objects:
##   - BamTallyParam
##   - VariantTallySummaryParam (defined in gmapR for now)
## - In theory, could make tallyVariants dispatch on both the 'x' and 'param'.
##   This would allow totally new tallying algorithms.

VariantTallyParam <- function(genome,
                              readlen = NA,
                              cycle_flank_width = 10L,
                              cycle_breaks = flankingCycleBreaks(readlen,
                                cycle_flank_width),
                              high_base_quality = 0L,
                              minimum_mapq = 13L, 
                              variant_strand = 1L, ignore_query_Ns = TRUE,
                              ...)
{
  if (variant_strand < 1 || variant_strand > 2)
    stop("'variant_strand' must be 1 or 2")
  if (!isTRUE(ignore_query_Ns))
    stop("'ignore_query_Ns' must be TRUE")
  args <- list(genome = genome,
               variant_strand = variant_strand,
               ignore_query_Ns = ignore_query_Ns,
               minimum_mapq = minimum_mapq,
               cycle_breaks = cycle_breaks,
               high_base_quality = high_base_quality,
               ...)
  do.call(BamTallyParam, args) 
}

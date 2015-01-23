### =========================================================================
### Tallying Alignments for Variant Calling
### -------------------------------------------------------------------------

### Thoughts on a tally abstraction:
###
### Want to support both bam_tally and applyPileups. Ideally, there
### would be one generic called tallyVariants(). We could dispatch on
### the parameter object, with a BamTallyVariantsParam() and a
### PileupTallyVariantsParam(). BamTallyVariantsParam could be
### provided by the gmapR package, which would depend on VariantTools,
### instead of the other way around.

setGeneric("tallyVariants", function(x, ...) standardGeneric("tallyVariants"))

defaultBPPARAM <- function() registered()[[1]]

setMethod("tallyVariants", "BamFile",
          function(x, param = TallyVariantsParam(...), ...,
                   BPPARAM = defaultBPPARAM())
          {
            if (!missing(param) && length(list(...)) > 0L) {
              warning("arguments in '...' are ignored when passing 'param'")
            }
            tally_region <- function(x, which, param) {
              iit <- bam_tally(x, param@bamTallyParam, which = which)
              ans <- variantSummary(iit, param@read_pos_breaks,
                                    param@high_base_quality,
                                    param@bamTallyParam@variant_strand == 0L,
                                    param@read_length)
              ## usage of start() is intential to avoid dropping indels
              ## that extend outside of window
              ans <- ans[start(ans) >= start(which) & start(ans) <= end(which)]
              if (!param@keep_extra_stats)
                mcols(ans) <- NULL
              ans[!(ans %over% param@mask)]
            }
            which <- param@bamTallyParam@which
            if (length(which) == 0L) {
              which <- tileGenome(seqlengths(param@bamTallyParam@genome),
                                  bpworkers(BPPARAM))
              which <- unlist(which, use.names=FALSE)
            }
            which <- as.list(which)
            ans <- bplapply(which, tally_region, x = x, param = param,
                            BPPARAM = BPPARAM)
            do.call(c, unname(ans))
          })

setMethod("tallyVariants", "BamFileList", function(x, ...) {
  stackSamples(VRangesList(bplapply(x, tallyVariants, ...)))
})

setMethod("tallyVariants", "character", function(x, ...) {
  tallyVariants(BamFile(x), ...)
})

setClass("TallyVariantsParam",
         representation(bamTallyParam = "BamTallyParam",
                        read_pos_breaks = "integerORNULL",
                        high_base_quality = "integer",
                        mask = "GenomicRanges",
                        keep_extra_stats = "logical",
                        read_length = "integer"))

TallyVariantsParam <- function(genome,
                               read_pos_breaks = NULL,
                               high_base_quality = 0L,
                               minimum_mapq = 13L,
                               variant_strand = 1L,
                               ignore_query_Ns = TRUE,
                               ignore_duplicates = TRUE,
                               mask = GRanges(),
                               keep_extra_stats = TRUE,
                               read_length = NA_integer_,
                               ...)
{
  if (!isSingleNumber(variant_strand) || !(variant_strand %in% c(0, 1, 2)))
    stop("'variant_strand' must be either 0, 1, or 2")
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
      mask = mask, keep_extra_stats = as.logical(keep_extra_stats),
      read_length = read_length)
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

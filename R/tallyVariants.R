### =========================================================================
### Tallying Alignments for Variant Calling
### -------------------------------------------------------------------------

setGeneric("tallyVariants", function(x, ...) standardGeneric("tallyVariants"))

defaultBPPARAM <- function() registered()[[1]]

setMethod("tallyVariants", "BamFile",
          function(x, param = TallyVariantsParam(...), ...,
                   BPPARAM = defaultBPPARAM())
          {
            tally_region <- function(x, ind, param) {
              which <- which[ind]
              iit <- bam_tally(x, param@bamTallyParam, which = which)
              ans <- variantSummary(iit, param@read_pos_breaks,
                                    param@high_base_quality,
                                    param@bamTallyParam@variant_strand == 0L)
              if (!param@keep_extra_stats)
                mcols(ans) <- NULL
              ans[!(ans %over% param@mask)]
            }
            which <- param@bamTallyParam@which
            if (length(which) == 0L) {
              which <- as(seqinfo(param@bamTallyParam@genome), "GenomicRanges")
            }
            ind <- seq_len(length(which))
            ans <- bplapply(ind, tally_region, x = x, param = param,
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
                        read_pos_breaks = "integer",
                        high_base_quality = "integer",
                        mask = "GenomicRanges",
                        keep_extra_stats = "logical"))

TallyVariantsParam <- function(genome,
                               read_pos_breaks = NULL,
                               high_base_quality = 0L,
                               minimum_mapq = 13L,
                               variant_strand = c(1L, 0L, 2L),
                               ignore_query_Ns = TRUE,
                               ignore_duplicates = TRUE,
                               mask = GRanges(),
                               keep_extra_stats = TRUE,
                               ...)
{
  variant_strand <- match.arg(variant_strand)
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
      mask = mask, keep_extra_stats = as.logical(keep_extra_stats))
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

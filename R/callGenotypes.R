### =========================================================================
### Genotype Calling
### -------------------------------------------------------------------------
###
### Routines for calling and manipulating genotypes.
###

setGeneric("callGenotypes",
           function(variants, cov, ...) standardGeneric("callGenotypes"))

### Strategy: calculate binomial probabilities and call
### whichever PL is 0 (p=1.0).

### This is equivalent to the LRT, as long as:
## p(0/0) is asking if the rate is the error rate (e)
## p(0/1) is asking if the rate is 0.5
## p(1/1) is asking if the rate is 1-e

## Generated fields:
## totalDepth: average coverage (DP)
## MIN_DP: minimum coverage
## GQ: genotype quality (p(GT) / p(AA) + p(RA) + p(AA))
## PL: likelihood of each genotype (RR, RA, AA) - binomial

phred <- function(p, max = 9999L) {
  ans <- pmin(round(-10 * log10(p)), max)
  storage.mode(ans) <- "integer"
  ans
}

compute_GQ <- function(gl, max = 99L) {
  phred(1-gl[cbind(seq_len(nrow(gl)), max.col(gl))] / rowSums(gl), max=max)
}

compute_GL <- function(ad, dp, p.error = 0.05) {
  gl <- cbind(dbinom(ad, dp, p.error),
              dbinom(ad, dp, 0.5),
              dbinom(ad, dp, 1-p.error))
  gl * 1.0 / gl[cbind(seq_len(nrow(gl)), max.col(gl, ties.method="first"))]
}

compute_GT <- function(gl) {
  c("0/0", "0/1", "1/1")[max.col(gl, ties.method="first")]
}

GenotypeRunVRanges <- function(ranges, proto, which, genome) {
  runs.gr <- GRanges(seqnames(which), shift(ranges, start(which) - 1L))
  ref <- getSeq(genome, resize(runs.gr, 1L))
  runs.vr <- VRanges(seqnames(runs.gr), ranges(runs.gr),
                     ref=ref, alt="<NON_REF>",
                     softFilterMatrix=matrix(NA, length(runs.gr),
                       ncol(softFilterMatrix(proto))),
                     seqlengths = seqlengths(proto))
  ml <- rep(list(rep(NA, length(runs.gr))), ncol(mcols(proto)))
  mcols(runs.vr) <- as(setNames(ml, names(mcols(proto))), "DataFrame")
  seqinfo(runs.vr) <- seqinfo(proto)
  runs.vr
}

## We take the GQ binning approach of GATK. An alternative strategy
## (from gvcftools) would be requiring data to fall into [x,y] where
## y <= x+max(3,x*0.3).

## We also take the average coverage for DP, like GATK, despite the
## gVCF spec recommending we take the min. Instead, the min is MIN_DP.

loadAndComputeGenotypesForRegion <- function(vcf, cov, which, ...) {
  vr <- readVcfAsVRanges(vcf, which=which)
  if (length(runValue(sampleNames(vr))) > 1L) {
    stop("variants must be from a single sample")
  }
  callGenotypesOneRegion(vr, cov, which, ...)
}

callGenotypesOneRegion <- function(variants, cov, which, param) {
  gl <- compute_GL(altDepth(variants), totalDepth(variants), param@p.error)
  variants$GT <- compute_GT(gl)
  variants$GQ <- compute_GQ(gl)
  variants$PL <- phred(gl)
  variants$MIN_DP <- rep(NA_integer_, length(variants))
  if (!is.null(cov)) {
    variants <- addGenotypeRuns(variants, cov, which, param)
  }
  variants
}

addGenotypeRuns <- function(variants, cov, which, param) {
  if (length(runLength(seqnames(variants))) > 1L) {
    stop("currently all ranges must be on the same sequence")
  }
  which <- keepSeqlevels(which, as.character(seqnames(which)))
  cov <- import(cov, which=which, as="NumericList")[[1L]]
  cov.v <- as.integer(cov)
  gl <- compute_GL(0, cov.v, param@p.error)
  gq <- compute_GQ(gl)
  gq.runs <- ranges(Rle(cut(gq, param@gq.breaks)))
  runs <- setdiff(gq.runs, shift(ranges(variants), 1L - start(which)))
  runs.vr <- GenotypeRunVRanges(runs, variants, which, param@genome)
  runs.vr$GT <- Rle("0/0", length(runs.vr))
  runs.vr$GQ <- viewMins(Views(gq, runs))
  runs.vr$PL <- matrix(apply(phred(gl), 2,
                             function(gli) viewMins(Views(gli, runs))),
                       ncol=3L)
  cov.views <- Views(cov.v, runs)
  totalDepth(runs.vr) <- viewMeans(cov.views)
  runs.vr$MIN_DP <- viewMins(cov.views)
  c(variants, runs.vr)
}

setClass("CallGenotypesParam",
         representation(gq.breaks = "numeric",
                        p.error = "numeric",
                        genome = "ANY",
                        which = "GenomicRangesList"))

CallGenotypesParam <- function(genome, gq.breaks = c(0, 5, 20, 60, Inf),
                               p.error = 0.05,
                               which = tileGenome(seqinfo(genome), ntile=ntile),
                               ntile = 100L)
{
  if (!hasMethod("getSeq", class(genome))) {
    if (!requireNamespace("gmapR")) {
        stop("no getSeq() method for 'genome' ",
             "and the gmapR package is not installed to convert it to ",
             "a GmapGenome object")
    }
    genome <- gmapR::GmapGenome(genome)
  }
  if (any(is.na(gq.breaks)) || any(gq.breaks < 0)) {
    stop("'gq.breaks' values must be non-negative and non-NA")
  }
  if (!isSingleNumber(p.error) || p.error < 0) {
    stop("'p.error' must be a single, non-negative, non-NA number")
  }
  if (!isSingleNumber(ntile) || ntile < 0) {
    stop("'ntile' must be a single, non-negative, non-NA number")
  }
  new("CallGenotypesParam", gq.breaks=gq.breaks, p.error=p.error,
      genome=genome, which=which)
}

setMethod("callGenotypes", c("VRanges", "BigWigFile"),
          function(variants, cov,
                   param = CallGenotypesParam(variants),
                   BPPARAM = defaultBPPARAM())
          {
            sample <- runValue(sampleNames(variants))
            if (length(sample) > 1L) {
              stop("variants must be from a single sample")
            }
            which <- unlist(param@which)
            f <- factor(findOverlaps(variants, which, select="arbitrary"),
                        seq_len(length(which)))
            vl <- split(variants, f)
            ansl <- bpmapply(callGenotypesOneRegion, vl, as.list(which),
                             MoreArgs=list(
                               cov=cov,
                               param=param
                               ),
                             BPPARAM=BPPARAM)
            ans <- do.call(c, unname(ansl))
            sampleNames(ans) <- sample
            ans
          })

setMethod("callGenotypes", c("TabixFile", "BigWigFile"),
          function(variants, cov,
                   param = CallGenotypesParam(variants),
                   BPPARAM = defaultBPPARAM())
          {
            which <- as.list(unlist(param@which))
            ansl <- bplapply(which, loadAndComputeGenotypesForRegion,
                             vcf=variants, bigwig=cov, param=param,
                             BPPARAM=BPPARAM)
            ans <- do.call(c, unname(ansl))
            sample <- na.omit(unique(sampleNames(ans)))
            sampleNames(ans) <- sample
            ans
          })

setMethod("show", "CallGenotypesParam", function(object) {
  cat("A", class(object), "object\n", sep = " ")
  cat(gmapR:::showSlots(object, count = FALSE), sep = "")
})

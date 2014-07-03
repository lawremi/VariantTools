### =========================================================================
### Genotype Calling
### -------------------------------------------------------------------------
###
### Routines for calling and manipulating genotypes.
### 

### Genotype representation
## Currently, VRanges supports the conventional VCF genotype encoding
## "X/X" in the "GT" metadata column. This seems OK for now. One issue
## is positions that are heterozygous for two alt alleles. VRanges has
## only a single alt per row, so in these cases each row would have
## "./1" as its genotype. After conversion to VCF, we could have a
## special collapse function that merges these, if needed.

### Genotype calling

## callGT <- function(x, min.het.freq = 0.2, min.hom.freq = 1 - min.het.freq,
##                    callable.cov = 20L)
## {
##   ifelse(altFraction(x) < min.het.freq,
##          ifelse(refFraction(x) > min.hom.freq,
##                 ifelse(totalDepth(x) >= callable.cov,
##                        "0/0",
##                        "0/."),
##                 ifelse(refFraction(x) >= min.het.freq,
##                        "0/o", # 0/o ('o' means other)
##                        "o/o") # o/o
##                 ),
##          ifelse(altFraction(x) > min.hom.freq,
##                 ifelse(totalDepth(x) >= callable.cov,
##                        "1/1",
##                        "./1"),
##                 ifelse(refFraction(x) >= min.het.freq,
##                        "0/1",
##                        "o/1"))) # o/1
## }

## addGenotypeRuns <- function(x, cov, callable.cov, genome)
## {
##   callable.runs <- slice(cov, callable.cov)
##   x.rl <- as(x, "RangesList")
##   wt.runs <- setdiff(callable.runs, x.rl)
##   nocall.runs <- setdiff(gaps(ranges(callable.runs)), x.rl)
##   runsToGRanges <- function(runs, gt) {
##     gr <- as(runs, "GRanges")
##     ref <- getSeq(genome, resize(gr, 1L))
##     vr <- VRanges(seqnames(gr), ranges(gr), ref=ref,
##                   totalDepth=viewMins(Views(cov, runs)),
##                   softFilterMatrix=matrix(NA, length(gr),
##                     ncol(softFilterMatrix(x))))
##     ml <- rep(list(rep(NA, length(gr))), ncol(mcols(x)))
##     mcols(vr) <- as(setNames(ml, names(mcols(x))), "DataFrame")
##     mcols(vr)$GT <- gt
##     vr
##   }
##   c(x,
##     runsToGRanges(wt.runs, "0/0"),
##     runsToGRanges(nocall.runs, "./."))
## }

setGeneric("callGenotypes", function(x, ...) standardGeneric("callGenotypes"))
## setMethod("callGenotypes", "VRanges",
##           function(x, cov = NULL, p.het = 0.5, p.error = 0.05, power = 0.99,
##                    genome = GmapGenome(genome(x)))
##           {
##             if (!is.numeric(p.het) || is.na(p.het)) {
##               stop("'p.het' must be a non-NA number")
##             }
##             if (p.het < 0 || p.het > 0) {
##               stop("'p.het' must be in [0,1]")
##             }
##             min.het.freq <- lrtFreqCutoff(p.het, p.error)
##             filters <- VariantCallingFilters(p.lower=p.het, p.error=p.error)
##             callable.cov <- minCallableCoverage(filters, power=power)            
##             x$GT <- callGT(x, min.het.freq, callable.cov)
##             if (!is.null(cov)) {
##               x <- addGenotypeRuns(x, cov, callable.cov, genome)
##             }
##             x
##           })

### Alternative strategy: calculate binomial probabilities and call
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

phred <- function(p) {
  round(ifelse(p == 1.0, 0, -10 * log10(p)))
}

compute_GQ <- function(gl, max = 99L) {
  pmin(phred(1-gl[cbind(seq_len(nrow(gl)), max.col(gl))] / rowSums(gl)), max)
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

GenotypeRunVRanges <- function(seqnames, ranges, proto, genome) {
  runs.gr <- GRanges(seqnames, ranges)
  ref <- getSeq(genome, resize(runs.gr, 1L))
  runs.vr <- VRanges(seqnames(runs.gr), ranges(runs.gr),
                     ref=ref, alt="<NON_REF>",
                     softFilterMatrix=matrix(NA, length(runs.gr),
                       ncol(softFilterMatrix(proto))))
  ml <- rep(list(rep(NA, length(runs.gr))), ncol(mcols(proto)))
  mcols(runs.vr) <- as(setNames(ml, names(mcols(proto))), "DataFrame")
  runs.vr
}

## We take the GQ binning approach of GATK. An alternative strategy
## (from gvcftools) would be requiring data to fall into [x,y] where
## y <= x+max(3,x*0.3).

## We also take the average coverage for DP, like GATK, despite the
## gVCF spec recommending we take the min. Instead, the min is MIN_DP.

addGenotypeRuns <- function(x, cov, gq.breaks, p.error, genome) {
  computeRunsForChr <- function(chr) {
    gl <- compute_GL(0, as.vector(cov[[chr]]), p.error)
    gq <- compute_GQ(gl)
    gq.runs <- ranges(Rle(cut(gq, gq.breaks)))
    runs <- setdiff(gq.runs, ranges(x)[seqnames(x) == chr])
    runs.vr <- GenotypeRunVRanges(chr, gq.runs, x, genome)
    runs.vr$GT <- Rle("0/0", length(runs.vr))
    runs.vr$GQ <- viewMins(Views(gq, runs))
    runs.vr$PL <- apply(phred(gl), 2, function(gli) viewMins(Views(gli, runs)))
    cov.views <- Views(cov[[chr]], ranges(gq.views))
    totalDepth(vr) <- viewMeans(cov.views)
    runs.vr$MIN_DP <- viewMins(cov.views)
    runs.vr
  }
  runs <- do.call(c, unname(lapply(names(cov), computeRunsForChr)))
  c(x, runs)
}

setMethod("callGenotypes", "VRanges",
          function(x, cov = NULL, gq.breaks = c(0, 5, 20, 60, Inf),
                   p.error = 0.05, genome = GmapGenome(unique(genome(x))))
                   {
                     gl <- compute_GL(altDepth(x), totalDepth(x), p.error)
                     x$GT <- compute_GT(gl)
                     x$GQ <- compute_GQ(gl)
                     x$PL <- phred(gl)
                     if (!is.null(cov)) {
                       x <- addGenotypeRuns(x, cov, gq.breaks, p.error, genome)
                     }
                     x
                   })

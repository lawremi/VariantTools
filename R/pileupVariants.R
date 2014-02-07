### =========================================================================
### Generating Pileups for Variant Calling
### -------------------------------------------------------------------------

## TODO: add calculation of mean base quality
## TODO: support the BPPARAM -- how? can we pass it to applyPileups?

pileupVariants <- function(bams, genome, param = PileupParam(),
                           minAltDepth = 1L, baseOnly = TRUE,
                           BPPARAM = defaultBPPARAM(), ...)
{
  pileupFiles <- PileupFiles(bams)
  plpWhat(param) <- union(plpWhat(param), "seq")
  bamNames <- names(bams)
  if (is.null(bamNames))
    bamNames <- seq_len(length(bams))
  pileupFun <- function(x) {
    tally <- x$seq
    sn <- Rle(names(x$seqnames), x$seqnames)
    pos.gr <- GRanges(sn, IRanges(x$pos, width=1))
    seq.gr <- pos.gr
    seqnameStyle(seq.gr) <- seqnameStyle(genome)
    ref.pos <- getSeq(genome, seq.gr, as.character=TRUE)
    ref <- rep(ref.pos, each=nrow(tally)*ncol(tally))
    bases <- rep(rownames(tally), ncol(tally) * dim(tally)[3])
    is.ref <- ref == bases
    nalt <- (nrow(tally)-1)
    per.pos <- nalt*ncol(tally)
    per.samp.pos <- rep(nalt, ncol(tally) * dim(tally)[3])
    sampleNames <- rep(Rle(bamNames, rep(nalt, length(bamNames))), dim(tally)[3])
    vr <- VRanges(rep(sn, each=per.pos),
                  IRanges(rep(x$pos, each=per.pos), width=1),
                  alt = bases[!is.ref],
                  ref = rep(ref[is.ref], each=per.pos),
                  altDepth = tally[!is.ref],
                  refDepth = Rle(tally[is.ref], per.samp.pos),
                  totalDepth = Rle(colSums(tally), per.samp.pos),
                  sampleNames = sampleNames)
    if (dropZeros) {
      vr <- vr[altDepth(vr) >= minAltDepth]
    }
    if (baseOnly) {
      vr <- vr[ref(vr) != "N" & alt(vr) != "N"]
    }
    vr
  }
  gr <- do.call(c, applyPileups(pileupFiles, pileupFun, param=param))
  gr
}

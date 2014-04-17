### =========================================================================
### Generating Pileups for Variant Calling
### -------------------------------------------------------------------------

## TODO: add calculation of mean base quality
## TODO: support the BPPARAM -- how? can we pass it to applyPileups?

pileupVariants <- function(bams, genome, param = ApplyPileupsParam(),
                           minAltDepth = 1L, baseOnly = TRUE,
                           BPPARAM = defaultBPPARAM())
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
    if (!all(seqlevels(seq.gr) %in% seqlevels(genome))) {
      seqlevelsStyle(seq.gr) <- seqlevelsStyle(genome)
    }
    ref.pos <- as.character(getSeq(genome, seq.gr))
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
                  ref = rep(ref[is.ref], each=nalt),
                  altDepth = tally[!is.ref],
                  refDepth = Rle(tally[is.ref], per.samp.pos),
                  totalDepth = Rle(colSums(tally), per.samp.pos),
                  sampleNames = sampleNames)
    vr <- vr[altDepth(vr) >= minAltDepth]
    if (baseOnly) {
      vr <- vr[ref(vr) != "N" & alt(vr) != "N"]
    }
    vr
  }
  gr <- do.call(c, applyPileups(pileupFiles, pileupFun, param=param))
  gr
}

### Some unexported work on a converter from a VRanges to a "pileup" GRanges:
pileupGRanges <- function(x) {
  indel <- isIndel(x)
  if (any(indel, na.rm=TRUE)) {
    x <- x[-which(indel)]
  }
  gr <- GRanges(seqnames(x), ranges(x), strand(x))
  sm <- selfmatch(gr)
  uniq <- sm == seq_len(length(gr))
  map <- integer()
  map[sm[uniq]] <- seq_len(sum(uniq))
  pos <- map[sm]
  base <- match(c(alt(x), ref(x)), DNA_BASES)
  depth <- c(altDepth(x), refDepth(x))
  samp <- as.factor(sampleNames(x))
  pileup <- array(0L, c(sum(uniq), length(levels(samp)), length(DNA_BASES)),
                  dimnames=list(NULL, levels(samp), base=DNA_BASES))
  pileup[cbind(pos, samp, base)[!is.na(base),]] <- depth[!is.na(base)]
  ugr <- gr[uniq]
  ugr$ref[pos] <- ref(x)
  ugr$pileup <- pileup
  ugr
}

pileupGRangesOverGenome <- function(x) {
  indel <- isIndel(x)
  if (any(indel, na.rm=TRUE)) {
    x <- x[-which(indel)]
  }
  genome <- unlist(tileGenome(seqinfo(x), tilewidth=1L))
  gr <- GRanges(seqnames(x), ranges(x), "*")
  m <- match(gr, genome)
  base <- match(c(alt(x), ref(x)), DNA_BASES)
  depth <- c(altDepth(x), refDepth(x))
  samp <- as.factor(sampleNames(x))
  pileup <- array(0L, c(length(genome), length(levels(samp)), length(DNA_BASES)),
                  dimnames=list(NULL, levels(samp), base=DNA_BASES))
  pileup[cbind(m, samp, base)[!is.na(base),]] <- depth[!is.na(base)]
  genome$REF <- "."
  genome$REF[m] <- ref(x)
  genome$REF <- DNAStringSet(genome$REF)
  genome$pileup <- pileup
  genome
}



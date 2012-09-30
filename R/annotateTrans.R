## Transcript-level variant summary:
## - callable fraction (see callableFraction)
## - functional variant count, essentially just tabulating tx ids
##

annotateTrans <- function(txdb, cov, anno_gr) {
  trans <- cdsBy(txdb, "tx", use.names=TRUE)
  values(trans)$callable.fraction <- callableFraction(txdb, cov, prob = 0.8)
  mcols(trans)$prot.alt.count <- table(factor(anno_gr$TXID, names(trans)))
  trans
}

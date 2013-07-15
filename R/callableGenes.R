setGeneric("callableFraction",
           function(tx, reads, ...) standardGeneric("callableFraction"))

setMethod("callableFraction", c("GRangesList", "ANY"),
          function(tx, reads, ...) {
            callable <- callCallable(reads, ...)
            seqlevels(tx) <- names(cov)
            tx_gr <- unlist(tx, use.names = FALSE)
            tx_ranges <- split(ranges(tx_gr), seqnames(tx_gr))
            exon_sums <- unsplit(viewSums(Views(callable, tx_ranges)),
                                 seqnames(tx_gr))
            tx_sums <- rowsum(exon_sums, names(tx)[togroup(tx)])
            tx_sums / sum(width(tx))
          })

setMethod("callableFraction", c("TranscriptDb", "ANY"),
          function(tx, reads, ...) {
            callGeneric(exonsBy(tx), reads, ...)
          })

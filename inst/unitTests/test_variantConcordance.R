test_checkVariantConcordance() {
  p53 <- gmapR:::exonsOnTP53Genome("TP53")
  bams <- LungCancerLines::LungCancerBamFiles()
  bam <- bams$H1993
  tally.param <- VariantTallyParam(gmapR::TP53Genome(), 
                                   readlen = 100L,
                                   high_base_quality = 23L,
                                   which = range(p53))
  called.variants <- callVariants(bam, tally.param)
  concord <- checkVariantConcordance(called.variants,
                                     called.variants)
  checkTrue(1L, concord$fraction)
}

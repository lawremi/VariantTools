test_calculateConcordanceMatrix_fileInput <- function() {
  p53 <- gmapR:::exonsOnTP53Genome("TP53")
  bams <- LungCancerLines::LungCancerBamFiles()
  tally.param <- VariantTallyParam(gmapR::TP53Genome(), 
                                   readlen = 100L,
                                   high_base_quality = 23L,
                                   which = range(p53))
  called.variants1 <- callVariants(bams$H1993, tally.param)
  called.variants2 <- callVariants(bams$H2073, tally.param)
  
  varFile1 <- tempfile()
  varFile2 <- tempfile()
  save(called.variants1, file=varFile1)
  save(called.variants2, file=varFile2)
  on.exit(unlink(varFile1), add=TRUE)
  on.exit(unlink(varFile2), add=TRUE)
  variantFiles <- c(varFile1, varFile2)

  mat <- calculateConcordanceMatrix(variantFiles)
  expected <- matrix(c(1, 1, 1, 1), ncol = length(variantFiles))
  rownames(expected) <- colnames(expected) <- variantFiles
  checkIdentical(expected, mat)
}

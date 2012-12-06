test_getConcordantCliques_twoConcord() {
  input <- matrix(c(1, 0.9, 0.9, 1), ncol=2)
  threshold <- 0.8 
  output <- .getConcordantCliques(input, threshold)
  expectedOutput <- list("1" = c("id1", "id2"))
  checkIdentical(expectedOutput, output)
}

test_getConcordantCliques_noConcord() {
  input <- matrix(c(1, 0.1, 0.1, 1), ncol=2)
  threshold <- 0.8 
  output <- .getConcordantCliques(input, threshold)
  expectedOutput <- list("1" = "id1", "2" = "id2")
  checkIdentical(expectedOutput, output)
}

test_getConcordantCliques_twoIdenticalPairs() {
  ##A is C and B is D
  input <- matrix(c(1, 0, 0.9, 0,
                    0, 1, 0, .95,
                    .9, 0, 1, 0,
                    0, .95, 0, 1),
                  ncol=4)
  threshold <- 0.85
  rownames(input) <- colnames(input) <- LETTERS[seq_len(ncol(input))]
  output <- .getConcordantCliques(input, threshold)
  expectedOutput <- list("1" = c("A", "C"), "2" = c("B", "D"))
  checkIdentical(expectedOutput, output)
}

test_idVerify_twoConcordant() {
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
  concordReport <- idVerify(variantFiles, threshold=0.8)
  expected <- rep("concordant", 2L)
  names(expected) <- variantFiles
  checkIdentical(expected, concordReport)
}

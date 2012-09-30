plotFreq <- function(sam){ 
  gr1 <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/TN/",
                        sam, "/",
                        sam,
                        ".sample_specific_variants_granges.RData",
                        sep = "")))

  gr2 <- get(load(paste("/gne/research/data/cgp/2.0/",
                        sam,
                        "/rna_seq_merged/RData/",
                        sam,
                        ".sample_specific_variants_granges.RData",
                        sep = "")))
  over <- findOverlaps(gr1, gr2)
  from_ex <- gr1[queryHits(over)]
  from_rna <- gr2[subjectHits(over)]
  rat1 <- values(from_ex)$count/values(from_ex)$count.total
  rat2 <- values(from_rna)$count/values(from_rna)$count.total
  plot(rat1, rat2, pch = 19, xlim = c(0,1), ylim = c(0,1),
       xlab = "Exome read frequency", ylab = "RNA read frequency",
       main = paste("Somatic mutations in", sam, sep = " "))
  abline(0,1)
}

plotFreq2 <- function(sam){
  gr1 <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/Exome/",
                        sam, "/results/",
                        sam,
                        ".filtered_variants_granges.RData",
                        sep = "")))
  gr2 <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/RNA/",
                        sam, "/results/",
                        sam,
                        ".variants_granges.RData",
                        sep = "")))
  SNPs <- get(load("~/svn/Packages/VariantTools/inst/extdata/exon_snps_ill25.RData"))
  over <- findOverlaps(gr1, SNPs)
  in_snp1 <- gr1[queryHits(over)]
  in_snp1 <- in_snp1[values(in_snp1)$count >5]
  over <- findOverlaps(gr2, SNPs)
  in_snp2 <- gr2[queryHits(over)]
  in_snp2 <- in_snp2[values(in_snp2)$count >5]

  over <- findOverlaps(in_snp1, in_snp2)
  from_ex <- in_snp1[queryHits(over)]
  from_rna <- in_snp2[subjectHits(over)]
  rat1 <- values(from_ex)$count/values(from_ex)$count.total
  rat2 <- values(from_rna)$count/values(from_rna)$count.total
  plot(rat1, rat2, pch = 19, xlim = c(0,1), ylim = c(0,1),
       xlab = "Exome read frequency", ylab = "RNA read frequency",
       main = paste("Exonic SNP positions in", sam, sep = " "), cex=0.1)
  abline(0,1)
} 

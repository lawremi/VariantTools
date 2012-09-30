plotFreq <- function(sam){
  par(mfrow = c(1,2))
  vep <- read.table(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/TN/",
                          sam,
                          "_variant_annotations.txt",sep = ""),
                    sep = "\t")
  rna_gr <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/RNA/",
                           sam,
                           "/results/",
                           sam,
                           ".variants_granges.RData", sep = "")
                     ))
  ex_gr <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/Exome/",
                          sam,
                          "/results/",
                          sam,
                          ".filtered_variants_granges.RData", sep = "")
                    ))
  rna_gr <- rna_gr[values(rna_gr)$count.total >10]
  ns <- vep[vep$V7 == "NON_SYNONYMOUS_CODING",]
  s <- vep[vep$V7 == "SYNONYMOUS_CODING",]
  s_locs <- paste("chr", s$V2, ":*", sep = "")
  ns_locs <- paste("chr", ns$V2, ":*", sep = "")
  ns_ex_gr <- ex_gr[values(ex_gr)$location %in% ns_locs]
  s_ex_gr <- ex_gr[values(ex_gr)$location %in% s_locs]
  ns_over <- findOverlaps(ns_ex_gr, rna_gr)
  s_over <- findOverlaps(s_ex_gr, rna_gr)
  ns_over_ex <- ns_ex_gr[queryHits(ns_over)]
  ns_over_rna <- rna_gr[subjectHits(ns_over)]
  s_over_ex <- s_ex_gr[queryHits(s_over)]
  s_over_rna <- rna_gr[subjectHits(s_over)]
  plot(values(ns_over_ex)$count/values(ns_over_ex)$count.total, values(ns_over_rna)$count/values(ns_over_rna)$count.total, pch = 19, cex = 0.4, col = "blue", xlim = c(0,1), ylim =c(0,1), ylab = "Alternate read frequency in RNA", xlab = "Alternate read frequency in Exome", main = sam)
  par(new=T)
  plot(values(s_over_ex)$count/values(s_over_ex)$count.total, values(s_over_rna)$count/values(s_over_rna)$count.total, pch = 19, cex = 0.4, col = "red", xlim = c(0,1), ylim = c(0,1), xlab= NA, ylab = NA)
  abline(0,1)
  v2 <- (values(ns_over_ex)$count/values(ns_over_ex)$count.total)/(values(ns_over_rna)$count/values(ns_over_rna)$count.total)
  v1 <- (values(s_over_ex)$count/values(s_over_ex)$count.total)/(values(s_over_rna)$count/values(s_over_rna)$count.total)
  breaks = seq(-5, 20, by = 0.5)
  names = breaks[1:(length(breaks)-1)]
  c1 = hist(log2(v1), breaks = breaks, plot = FALSE)$count
  c2 = hist(log2(v2), breaks = breaks, plot = FALSE)$count
  barplot(rbind(c1/sum(c1), c2/sum(c2)), beside =T, legend = c("SYN", "Non-SYN"), names = names, las=2, ylab = "Normalized fraction", xlab = "log2 (exome freq / rna freq)" )
}




getGenes <- function(sam, freq){
  vep <- read.table(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/TN/",
                          sam,
                          "_variant_annotations.txt",sep = ""),
                    sep = "\t")
  if(file.exists(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/RNA/",
                       sam,
                       "/results/",
                       sam,
                       ".variants_granges.RData", sep = "")
                 )){
    rna_gr <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/RNA/",
                             sam,
                             "/results/",
                             sam,
                             ".variants_granges.RData", sep = "")
                       ))
    ex_gr <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/Exome/",
                            sam,
                            "/results/",
                            sam,
                            ".filtered_variants_granges.RData", sep = "")
                      ))
  }else{
    return(NULL)
  }
  rna_gr <- rna_gr[values(rna_gr)$count.total >10]
  ns <- vep[vep$V7 == "NON_SYNONYMOUS_CODING",]
  s <- vep[vep$V7 == "SYNONYMOUS_CODING",]
  s_locs <- paste("chr", s$V2, ":*", sep = "")
  ns_locs <- paste("chr", ns$V2, ":*", sep = "")
  ns_ex_gr <- ex_gr[values(ex_gr)$location %in% ns_locs]
  s_ex_gr <- ex_gr[values(ex_gr)$location %in% s_locs]
  ns_over <- findOverlaps(ns_ex_gr, rna_gr)
  s_over <- findOverlaps(s_ex_gr, rna_gr)
  ns_over_ex <- ns_ex_gr[queryHits(ns_over)]
  ns_over_rna <- rna_gr[subjectHits(ns_over)]
  s_over_ex <- s_ex_gr[queryHits(s_over)]
  s_over_rna <- rna_gr[subjectHits(s_over)]
  ns_over_rna[values(ns_over_rna)$count/values(ns_over_rna)$count.total >freq]->high
  vep[which(ns_locs %in% values(high)$location),]->damaging
  split <-strsplit(as.character(damaging$V14), split = "=")
  gene_ids2 <- lapply(split, function(x) x[5])
  gene_ids2 <- unique(unlist(gene_ids2))
  return(gene_ids2)
}

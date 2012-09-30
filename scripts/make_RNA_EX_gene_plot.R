library(VariantTools)
library(CGPtools)
library(parallel)
library(TxDb.Hsapiens.BioMart.igis)
library(org.Hs.eg.db)


txdb <- TxDb.Hsapiens.BioMart.igis
genes <- cdsBy(txdb, "gene")

x <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])


##get_sams
meta<-getCGPMetadata(JSON="http://research.gene.com/cgp/release_2.1/released.json")
colon <- c('Colon', 'Colon Sigmoid', 'Rectum', 'Cecum', 'Ileum', 'Intestine Small')
exome <- meta[meta$project == "exome",]
normals <- as.character(exome$srcid[exome$tissue_diagnosis %in% c("Normal", "Normal Adjacent Tumor")])
tissues <- as.character(exome$primary_tissue)
names(tissues) <- as.character(exome$srcid)
tissues[normals]<-"normal"
tissues[tissues %in% colon]<-"Colon"
colon_sams <- names(tissues[tissues=="Colon"])
colon_sams <- as.list(colon_sams)
##
gene <- "TP53"
gene_id <- unlist(mget(gene, org.Hs.egSYMBOL2EG))
gene_gr <- genes[[paste("GeneID:", gene_id,sep = "")]]

list <- mclapply(colon_sams, mc.cores=detectCores(), function(x){freq_info(x, gene_gr)})
x <- do.call(c, list)
pdf(file = "TP53_EXvR_freq.pdf", height=8, width=8)
plot(values(x)$count/values(x)$count.total, values(x)$norm_alt/values(x)$norm_cov,
     xlim = c(0,1), ylim = c(0,1), pch = 19, xlab = "Alt read frequency Exome",
     ylab = "Alt read frequency RNA", main =gene )
abline(0,1)
dev.off()
freq_info <- function(sam, gene_gr){
  f1 <- paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/TN/",
              sam,
              "/",
              sam,
              ".sample_specific_variants_granges.RData", sep = "")
  f2 <- paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/RNA/",
              sam,
              "/results/",
              sam,
              ".variants_granges.RData", sep = "")
  f3 <- paste("/gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS115/",
              sam,
              "/results/",
              sam
              ,".coverage.RData", sep = "")
  if(file.exists(f1) & file.exists(f2) &file.exists(f3)){
    gr <- get(load(f1))
    raw<- get(load(f2))
    cov <- get(load(f3))
    gr_in_gene <- findOverlaps(gr, gene_gr)
    gr <- gr[queryHits(gr_in_gene)]
    if(length(gr)==0){
      return(GRanges())
    }
    gr<- sort(gr)
    raw<-sort(raw)
    values(gr)$norm_ref <- as.numeric(NA)
    values(gr)$norm_alt <- as.numeric(0)
    over2 <- which(paste(values(gr)$location,values(gr)$read) %in% paste(values(raw)$location, values(raw)$read))
    over <- which(paste(values(raw)$location,values(raw)$read) %in% paste(values(gr)$location, values(gr)$read))
    if(length(over)==0 |length(over2) ==0){
      return(GRanges())
    }
    values(gr[over2])$norm_ref <- values(raw[over])$count.ref
    values(gr[over2])$norm_alt <- values(raw[over])$count
    ##    values(gr[-over2])$norm_alt <- 0
    
    grl <- split(gr, as.factor(seqnames(gr)))
    names <- as.list(names(grl))
    n <- lapply(names, function(x) { as.numeric(cov[[x]][ranges(grl[[x]])]) })
    n <- unlist(n)
    gr <- unlist(grl)
    values(gr)$norm_cov<-n 
    return(gr)
  }else{
    return(GRanges())   
  }
}

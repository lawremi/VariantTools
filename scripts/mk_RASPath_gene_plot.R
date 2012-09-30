library(CGPtools)
library(VariantTools)
library(parallel)
library(TxDb.Hsapiens.BioMart.igis)
library(org.Hs.eg.db)
source("~/svn/Packages/VariantTools/scripts/getAlleleByLocus_mod.R")
txdb <- TxDb.Hsapiens.BioMart.igis
genes <- cdsBy(txdb, "gene")

x <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

meta<-getCGPMetadata(JSON="http://research.gene.com/cgp/release_2.1/released.json")

data <- get(load("~/svn/Projects/CRCanalysis/trunk/data/RNAseq_64_matchedDFs_PTIDdf.rda"))
colon_sams <- as.list(as.character(data[,2]))


##APC
kp <- unlist(mget(c("APC"), org.Hs.egSYMBOL2EG))
apc <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic") 
kp <- unlist(mget(c("WNT5A"), org.Hs.egSYMBOL2EG))
wnt5a <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("CTNNB1"), org.Hs.egSYMBOL2EG))
beta <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("FZD1"), org.Hs.egSYMBOL2EG))
fzd1 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("FZD2"), org.Hs.egSYMBOL2EG))
fzd2 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("DVL1"), org.Hs.egSYMBOL2EG))
dvl1 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("DVL2"), org.Hs.egSYMBOL2EG))
dvl2 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("DVL3"), org.Hs.egSYMBOL2EG))
dvl3 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
apc_freq<-freq(apc)
wnt5a_freq<-freq(wnt5a)
beta_freq<-freq(beta)
fzd1_freq<-freq(fzd1)
fzd2_freq<-freq(fzd2)
dvl1_freq<-freq(dvl1)
dvl2_freq<-freq(dvl2)
dvl3_freq<-freq(dvl3) 
mat <- cbind(apc_freq, wnt5a_freq, beta_freq,fzd1_freq, fzd2_freq, dvl1_freq, dvl2_freq, dvl3_freq)
colname <- c("APC", "WNT5a", "beta-catenin", "Fzd1", "Fzd2", "dvl1", "dvl2", "dvl3")
mat2 <- mat[sort.list(mat[,1]), ]

par(mar= c(8,5, 3, 5))
image(mat2, col = tim.colors(64), axes = F, main = "WNT pathway mutations")
axis(1, at=seq(0,1, by=1/(length(rownames(mat))-1)),lab=rownames(mat2), las = 2)
axis(2, at=seq(0,1, by=1/(length(colname)-1)),lab=colname, las = 2)
image.plot(mat2, legend.only=T, legend.lab="Mutation Frequency", legend.mar=4.0)

kp <- unlist(mget(c("FAT4"), org.Hs.egSYMBOL2EG))
fat4 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("PIK3CA"), org.Hs.egSYMBOL2EG))
pik3ca <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("CSMD1"), org.Hs.egSYMBOL2EG))
csmd1 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("CSMD2"), org.Hs.egSYMBOL2EG))
csmd2 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("CSMD3"), org.Hs.egSYMBOL2EG))
csmd3 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("SYNE1"), org.Hs.egSYMBOL2EG))
syne1 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic") 
kp <- unlist(mget(c("MUC6"), org.Hs.egSYMBOL2EG))
muc6 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("ATM"), org.Hs.egSYMBOL2EG))
atm <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("OBSCN"), org.Hs.egSYMBOL2EG))
obscn <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
##SMAD4
kp <- unlist(mget(c("SMAD4"), org.Hs.egSYMBOL2EG))
smad4 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("SMAD3"), org.Hs.egSYMBOL2EG))
smad3 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("SMAD2"), org.Hs.egSYMBOL2EG))
smad2 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
##APC
kp <- unlist(mget(c("APC"), org.Hs.egSYMBOL2EG))
apc <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
##TP53
kp <- unlist(mget(c("TP53"), org.Hs.egSYMBOL2EG))
tp53 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic") 
##RAS
kp <- unlist(mget(c("KRAS"), org.Hs.egSYMBOL2EG))
kras <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("NRAS"), org.Hs.egSYMBOL2EG))
nras <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("HRAS"), org.Hs.egSYMBOL2EG))
hras <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
##RAF
kp <- unlist(mget(c("BRAF"), org.Hs.egSYMBOL2EG))
braf <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic") 
kp <- unlist(mget(c("RAF1"), org.Hs.egSYMBOL2EG))
craf <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic") 
kp <- unlist(mget(c("ARAF"), org.Hs.egSYMBOL2EG))
araf <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic") 
##MEK
kp <- unlist(mget(c("MAP2K1"), org.Hs.egSYMBOL2EG))
map2k1 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("MAP2K2"), org.Hs.egSYMBOL2EG))
map2k2 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
##ERK
kp <- unlist(mget(c("MAPK1"), org.Hs.egSYMBOL2EG))
mapk1 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic")
kp <- unlist(mget(c("MAPK3"), org.Hs.egSYMBOL2EG))
mapk3 <- getAlleleByLocus_mod(gene=kp, meta=meta, gene_db=genes, ids=colon_sams, tech="exome", type = "somatic") 

kras_freq <-freq(kras)
nras_freq <-freq(nras)
hras_freq <-freq(hras)
braf_freq <-freq(braf)
craf_freq <- freq(craf)
araf_freq<-freq(araf)
map2k1_freq <-freq(map2k1)
map2k2_freq <-freq(map2k2)
mapk1_freq <-freq(mapk1)
mapk3_freq <-freq(mapk3)
apc_freq<-freq(apc)
tp53_freq<-freq(tp53)
smad4_freq<-freq(smad4)
smad3_freq<-freq(smad3)
smad2_freq<-freq(smad2)
atm_freq<-freq(atm)
obscn_freq<-freq(obscn)
muc6_freq<-freq(muc6)
syne1_freq<-freq(syne1)
csmd1_freq<-freq(csmd1)
csmd2_freq<-freq(csmd2)
csmd3_freq<-freq(csmd3)
pik3ca_freq<-freq(pik3ca)
fat4_freq<-freq(fat4)

mat <- cbind(kras_freq, nras_freq, hras_freq, braf_freq, craf_freq,
             araf_freq, map2k1_freq, map2k2_freq, mapk1_freq,
             mapk3_freq, apc_freq, tp53_freq, smad4_freq, smad3_freq,
             smad2_freq, atm_freq, obscn_freq, syne1_freq, csmd1_freq,
             csmd2_freq, csmd3_freq,muc6_freq, pik3ca_freq, fat4_freq)

loh <- get(load("~/svn/Projects/CRCanalysis/trunk/data/loh_data.RData"))
cnv <- get(load("~/svn/Projects/CRCanalysis/trunk/data/cnv_data.RData"))
vec <- matrix(NA, nrow=64, ncol=2)
rownames(vec) <- rownames(mat)
colnames(cnv) %in% rownames(mat)->cols
cnv[paste("GeneID.", kp, sep = ""),cols]->kras_cn
loh[paste("GeneID.", kp, sep = ""),cols]->kras_loh
vec[names(kras_cn), 1]<-kras_cn
vec[names(kras_loh), 2]<-kras_loh

colname <- c("KRAS", "NRAS", "HRAS", "BRAF", "c-RAF", "a-RAF",
             "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "APC", "TP53",
             "SMAD4", "SMAD3", "SMAD2", "ATM", "OBSCN", "SYNE1",
             "CSMD1", "CSMD2", "CSMD3", "MUC6", "PIK3CA", "FAT4")

mat2 <- mat[sort.list(mat[,19]), ]
vec2 <- vec[sort.list(mat[,1]),]

pdf(file = "RAS_PATH_MUT.pdf", height = 8, width = 16)
#par(mfrow=c(2,1))
par(mar= c(8,5, 3, 5))
image(mat2, col = tim.colors(64), axes = F, main = "RAS pathway mutations")
axis(1, at=seq(0,1, by=1/(length(rownames(mat))-1)),lab=rownames(mat2), las = 2)
axis(2, at=seq(0,1, by=1/(length(colname)-1)),lab=colname, las = 2)
image.plot(mat2, legend.only=T, legend.lab="Mutation Frequency", legend.mar=4.0) 
dev.off()

pdf(file = "KRAS_LOH.pdf", height = 5, width = 16)
par(mar= c(8,5, 3, 5))
image(vec2, col = tim.colors(64), axes = F, main = "KRAS CN/LOH")
axis(1, at=seq(0,1, by=1/(64-1)),lab=rownames(mat2), las = 2)
axis(2, at=c(0.05, 0.95),lab=c("CNV state", "LOH state"), las = 2)
image.plot(vec2, legend.only=T, legend.lab="CN/LOH value", legend.mar=4.0)
dev.off() 

freq<- function(list){
  freq <- lapply(list, function(x) {
    if(length(x)>0){
      return(max(values(x)$count/values(x)$count.total))
    }else{
      return(0)
    }
  })
  return(unlist(freq))
}

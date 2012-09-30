library(CGPtools)
library(TxDb.Hsapiens.BioMart.igis)
library(parallel)
library(VariantTools)
source("~/svn/Packages/VariantTools/scripts/functional_vep.R")
txdb <- TxDb.Hsapiens.BioMart.igis 

colnames <- c("Uploaded_variation",
              "Location",
              "Allele", "Gene", "Feature",
              "Feature_type", "Consequence", "DNA_position",
              "CDS_position", "Protein_position",
              "Amino_acids", "Codons", "Existing_variation", "Extra")

dir <-"/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/TN/varEffOut/"
files <- dir(dir)
files <- files[grepl(pattern=".txt", files)]
file_sam <- sub("_variant_annotations.txt", "", files)

cov_dir <- "/gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS121/"


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

run_anno <- function(sam){
  file <- files[which(file_sam ==sam)]
  if(length(file) != 0){
    tab <- read.table(file.path(dir, file), sep = "\t", stringsAsFactors=FALSE, colClasses="character")
    colnames(tab)<-colnames
  }else{
    return(NULL)
  }
  vep_gr <- VEP2gr(tab)
  vep_gr <- getFunctional(vep_gr)
  cov_file <- paste(cov_dir, sam, "/results/", sam, ".coverage.RData", sep = "")
  if(file.exists(cov_file)){
    cov <- get(load(cov_file))
  }else{
    return(NULL)
  }
  anno <- annotateTrans(txdb, cov, vep_gr)
  return(anno)
}

all_anno <- lapply(colon_sams, function(sam){run_anno(sam)})
names <- names(all_anno[[1]])
names(all_anno) <- unlist(colon_sams)
mat <- matrix(ncol=73, nrow = 42969, NA)
colnames(mat) <- unlist(colon_sams)
rownames(mat) <- names
for(i in 1:73){
  if(!is.null(all_anno[[i]])){
    x <- values(all_anno[[i]])
    v1 <-ifelse(x$prot_alt_count>0, yes=1, no=NA)
    mat[,i] <- ifelse(x$callable_fraction>.8 & x$prot_alt_count==0, yes = 0, no=v1)
  }else{
  }
}

library(CGPtools)

colnames <- c("Uploaded_variation",
              "Location",
              "Allele", "Gene", "Feature",
              "Feature_type", "Consequence", "DNA_position",
              "CDS_position", "Protein_position",
              "Amino_acids", "Codons", "Existing_variation", "Extra")

dir <-"/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/TN/"
files <- dir(dir)
files <- files[grepl(pattern=".txt", files)]
file_sam <- sub("_variant_annotations.txt", "", files)

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


sam <- colon_sams[[1]]
DG <- damaging_genes(sam)

DG <- lapply(colon_sams, function(sam) damaging_genes(sam))

lung_DG <- lapply(lung_sams, function(sam) damaging_genes(sam))

SpG <- lapply(colon_sams, function(sam) splice_genes(sam))

damaging_genes <- function(SAM){
  x <- files[which(file_sam ==SAM)]
  if(length(x) != 0){
    tab <- read.table(file.path(dir, x), sep = "\t", stringsAsFactors=FALSE, colClasses="character")
    colnames(tab)<-colnames
    damaging <- tab[grepl(pattern="STOP_GAINED", tab$Consequence),]
    split <-strsplit(damaging$Extra, split = "=")
    gene_ids1 <- lapply(split, function(x) x[2])
    gene_ids1 <- unique(unlist(gene_ids1))
    damaging <- tab[(grepl(pattern="damaging", tab$Extra) | grepl(pattern="deleterious", tab$Extra)),]
    split <-strsplit(damaging$Extra, split = "=")
    gene_ids2 <- lapply(split, function(x) x[5])
    gene_ids2 <- unique(unlist(gene_ids2))
    gene_ids<-c(gene_ids1, gene_ids2)
    return(gene_ids)
  }else{
    return(NULL)
  }
}

splice_genes <- function(SAM){
  x <- files[which(file_sam ==SAM)]
  if(length(x) != 0){
    tab <- read.table(file.path(dir, x), sep = "\t", stringsAsFactors=FALSE, colClasses="character")
    colnames(tab)<-colnames
    damaging <- tab[grepl(pattern="SPLICE_SITE", tab$Consequence),]
    split <-strsplit(damaging$Extra, split = "=")
    gene_ids <- lapply(split, function(x) x[5])
    gene_ids <- unique(unlist(gene_ids))
    return(gene_ids)
  }else{
    return(NULL)
  }
}



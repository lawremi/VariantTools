library(genoset)
library(BigMatrix)
library(CGPtools)
library(VariantTools)
source("~/svn/Packages/VariantTools/scripts/getAlleleByLocus_mod.R")

meta<-getCGPMetadata(JSON="http://research.gene.com/cgp/released.json")
## ds = readGenoSet("/gne/research/data/cgp/SNPArray.V4/illumina_2.5M.genoset.genotypes.and.raw.RData")
ds = readCGPSNPArray(version="Genotype")   # Very recent addition
my.levels = levels(ds[ , , "geno"])

het.levels = c("AC","AG","AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG") # or something

locs_gr <- as(locData(ds), "GRanges")
seqlevels(locs_gr) <- paste("chr", seqlevels(locs_gr), sep = "")


sample.44.hets = ds[ , 44, "geno"] %in% het.levels
sample.56.hets = ds[ , 105, "geno"] %in% het.levels

###looking at positions in sample with sufficient coverage

FN <- mclapply(sams, mc.cores =detectCores(), function(sam) check_het(sam))

check_het <- function(SAM, return_gr =FALSE){
  cov <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS121/", SAM, "/results/",SAM,".coverage.RData", sep = "")))
  
  gr1 <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/",
                        SAM,
                        "/results/",
                        SAM,
                        ".filtered_variants_granges.RData", sep = "")))
  seqlevels(gr1) <- paste("chr", seqlevels(gr1), sep = "")
  names(cov) <- paste("chr", names(cov), sep = "")
  names <- seqlevels(locs_gr)
  grl <- split(locs_gr, as.factor(seqnames(locs_gr)))
  n <- lapply(names, function(x) { as.numeric(cov[[x]][ranges(grl[[x]])]) })
  values(locs_gr)$cov <- unlist(n)
  int <- which(pData(ds)$SAMPLE_ID ==sub("SAM", "", SAM))
  sample.hets = ds[ , int, "geno"] %in% het.levels
  het_locs <- locs_gr[sample.hets]
  high_cov_locs <- het_locs[values(het_locs)$cov >15]
  over <- findOverlaps(gr1, high_cov_locs)
  not_found <- high_cov_locs[-subjectHits(over)]
  ret <- 1- (length(not_found)/length(high_cov_locs))
  if(return_gr){
    gr2 <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/",
                          SAM,
                          "/results/",
                          SAM,
                          ".variants_granges.RData", sep = "")))
    seqlevels(gr2) <- paste("chr", seqlevels(gr2), sep = "")
   
    cat(ret," ", length(high_cov_locs),"\n", sep = "")
    over <- findOverlaps(gr2, not_found)
    gr3 <- gr2[queryHits(over)]
    var <- variantFilter(gr3,pval=0.001, useQual=TRUE, read_pos= "cycleCount.15.60")
    return(c(ret, length(high_cov_locs),
             length(not_found),
             length(var$raw_granges),
             var$low_count_rejects,
             var$frequency_rejects,
             var$fisher_test_rejects_2by2,
             var$read_pos_rejects))
  }else{
    return(c(ret, length(high_cov_locs)))
  }
}  

check <-getAlleleByLocus_mod(locus = not_found, meta=meta, ids=SAM, tech="exome", type = "raw")

var <- check[[1]][as.character(values(check[[1]])$ref) != as.character(values(check[[1]])$read)]
variantFilter(var)

gr2 <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/",
                      SAM,
                      "/results/",
                      SAM,
                      ".variants_granges.RData", sep = "")))    


** Actually it looks like the levels fell off of the PRD genotype
object.  Lovely.  They are

paste( c(rep("A",4),rep("C",4),rep("G",4),rep("T",4)),

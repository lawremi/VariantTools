plotT <- function(normal_sam, tumor_sam){
  dir <- "/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/Exome/"
  tumor_gr <- get(load(paste(dir, tumor_sam, "/results/", tumor_sam, ".filtered_variants_granges.RData", sep = "")))
  normal_gr <- get(load(paste(dir, normal_sam,"/results/", normal_sam, ".filtered_variants_granges.RData", sep = "")))
  normal_raw <- get(load(paste(dir, normal_sam, "/results/", normal_sam, ".variants_granges.RData", sep="")))
  normal_cov <-  get(load(paste("/gne/research/data/cgp/2.0/", normal_sam, "/exome_merged/RData/", normal_sam, ".cov.RData", sep = "")))
  VariantTools:::plotVarGenome(tumor_gr= tumor_gr, normal_gr = normal_gr, normal_raw=normal_raw, normal_cov=normal_cov)
}

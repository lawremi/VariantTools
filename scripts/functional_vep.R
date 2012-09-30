getFunctional<- function(vep_gr){
  vep_gr <- vep_gr[which(grepl(pattern = "GeneID", values(vep_gr)$Gene))]
  ns_inds <- which(grepl(pattern="NON_SYNONYMOUS_CODING",  values(vep_gr)$Consequence))
  stop_inds <- which(grepl(pattern="STOP_GAINED",  values(vep_gr)$Consequence))
  inds <- c(ns_inds, stop_inds)
  vep_gr <- vep_gr[inds]
  return(vep_gr)
}

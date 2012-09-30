RVW <- function(sam){
  f1 <- paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/TN/",
              sam,"/",
              sam,
              ".sample_specific_variants_granges.RData", sep = "")
  if(file.exists(f1)){
    gr <- get(load(f1))
  }else{
    return(NULL)
  }
  
  f2 <- paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/RNA/",
              sam,"/results/",
              sam, ".variants_granges.RData", sep = "")
  
  f3 <- paste("/gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS115/",
              sam, "/results/",
              sam, ".coverage.RData", sep = "")
  if(file.exists(f2)){
    
    rna <- get(load(f2))
    cov <- get(load(f3))
  }else{
    return(NULL)
  }
  val <- rnaValRate(gr=gr, rna=rna, cov=cov)
  return(val)
}

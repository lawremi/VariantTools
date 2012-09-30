out_dir <- "/gne/research/workspace/degenhj2/run_cgp3_ex/"
dir <- "/gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS121/"
sams <- dir(dir)

header <- "#!/bin/bash\nsource /gne/research/apps/modules/common/bashrc\nmodule load apps/R\n"
run <- "Rscript ~/svn/VariantTools/scripts/run_callVar_wrapper.R"
genome_info <-"hg19_IGIS21 /gne/research/apps/ngs_pipeline/dev/x86_64-linux-2.6-sles11/share 75"
out <- "/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/"
for(i in 1:length(sams)){
  if(file.exists(paste(dir, sams[i], "/bams/", sams[i], ".analyzed.bam", sep = ""))){
    bam_file <- paste(dir, sams[i], "/bams/", sams[i], ".analyzed.bam", sep = "")
    con <-file.path(out_dir, paste(sams[i], "_run.sh", sep = ""))
    writeLines(header, con = con)
    cat(paste(run, bam_file, genome_info, paste(out, sams[i], sep = ""), sams[i], "cgp_rerun\n"), append=T, file= con)
  }
  
}

out_dir <- "/gne/research/workspace/degenhj2/run_cgp3_rna/"
dir <- "/gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS115/"
sams <- dir(dir)

header <- "#!/bin/bash\nsource /gne/research/apps/modules/common/bashrc\nmodule load apps/R\n"
run <- "Rscript ~/svn/VariantTools/scripts/run_callVar_wrapper.R"
genome_info <-"hg19_IGIS21 /gne/research/apps/ngs_pipeline/dev/x86_64-linux-2.6-sles11/share 75"
out <- "/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/RNA/"
for(i in 1:length(sams)){
  if(file.exists(paste(dir, sams[i], "/bams/", sams[i], ".analyzed.bam", sep = ""))){
    bam_file <- paste(dir, sams[i], "/bams/", sams[i], ".analyzed.bam", sep = "")
    con <-file.path(out_dir, paste(sams[i], "_run.sh", sep = ""))
    writeLines(header, con = con)
    cat(paste(run, bam_file, genome_info, paste(out, sams[i], sep = ""), sams[i], "cgp_rerun\n"), append=T, file= con)
  }
  
}


library(CGPtools)
meta<-getCGPMetadata(JSON="http://research.gene.com/cgp/released.json")
meta <- meta[which(!is.na(meta$patient_num)),]
exome <- meta[meta$project=="exome",]
pat_ids <- unique(exome$patient_num)

header <- "#!/bin/bash\nsource /gne/research/apps/modules/common/bashrc\nmodule load apps/ngs_pipeline/3.4\nR_LIBS=/gne/home/degenhj2/R/x86_64-unknown-linux-gnu-library/2.14/:$R_LIBS\n"
run <- "Rscript ~/svn/VariantTools/scripts/run_TN_wrapper.R"
path <- "/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/Exome/"
out_dir <- "/gne/research/workspace/degenhj2/run_cgp2_ex_TN/"
out <- "/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/TN/"
for(i in 1:length(pat_ids)){
  ids <- as.character(exome[exome$patient_num ==pat_ids[i],"srcid"])

  tgr <- paste(path, ids[1], "/results/", ids[1], ".filtered_variants_granges.RData", sep = "")
  traw <- paste(path, ids[1], "/results/", ids[1], ".variants_granges.RData", sep = "")
  ngr <- paste(path, ids[2], "/results/", ids[2], ".filtered_variants_granges.RData", sep = "")
  nraw <- paste(path, ids[2], "/results/", ids[2], ".variants_granges.RData", sep = "")
  tcov <- paste("/gne/research/data/cgp/2.0/", ids[1], "/exome_merged/RData/", ids[1], ".cov.RData", sep = "")
  ncov <- paste("/gne/research/data/cgp/2.0/", ids[2], "/exome_merged/RData/", ids[2], ".cov.RData", sep = "")
  con1 <-file.path(out_dir, paste(ids[1], "_run_TN.sh", sep = ""))
  con2 <-file.path(out_dir, paste(ids[2], "_run_TN.sh", sep = ""))

  writeLines(header, con = con1)
  cat(paste(run, tgr, ngr, nraw, ncov,paste(out, ids[1], sep = ""), ids[1], "cgp_2.0_TN_rerun\n"), file=con1, append=TRUE)

  writeLines(header, con = con2)
  cat(paste(run, ngr, tgr, traw, tcov,paste(out, ids[2], sep = ""), ids[2], "cgp_2.0_TN_rerun\n"), file=con2, append=TRUE)
}


##########################################
## Running VariantEffectPredictor
#########################################

out_dir <- "/gne/research/workspace/degenhj2/run_cgp3_TN_vareff"
dir <- "/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/TN/"
sams <- dir(dir)
header <- "#!/bin/bash\nsource /gne/research/apps/modules/common/bashrc\nmodule load apps/igis\n"

for(i in 1:length(sams)){
  con = file.path(out_dir, paste(sams[i], "_varEff.sh", sep = ""))
  writeLines(header, con = con)
  cat(paste("zcat ",
            dir, "/",
            sams[i],
            "/",
            sams[i],
            "_sample_specific_variants.vcf.gz |variant_effect_predictor.pl --sift=b --polyphen=b --condel=b --hgnc --check_existing --output_file ",
            paste(dir,"varEffOut/", sep = ""),
            sams[i],
            "_variant_annotations.txt",
            sep = ""),
      file=con, append=TRUE)
              
 }

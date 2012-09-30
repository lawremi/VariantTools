library(VariantTools)
##library(RNASeqGenie)

cmd_args = commandArgs(trailingOnly=TRUE);
tumor_gr <- cmd_args[1];
normal_gr <- cmd_args[2];
normal_raw <- cmd_args[3];
normal_cov <- cmd_args[4];
save_dir <- cmd_args[5];
SAMid <- cmd_args[6];
project<-cmd_args[7];

normal_gr <- load(normal_gr)
normal_gr<-get(normal_gr)
tumor_gr<-load(tumor_gr)
tumor_gr<-get(tumor_gr)
normal_cov<- load(normal_cov)
normal_cov <- get(normal_cov)

if(file.exists(save_dir)){
  stop("The directory already exists")
}else{
  dir.create(dirname(save_dir))
  dir.create(save_dir)
}

TS<- tumorNormalCompare(tumor_gr=tumor_gr,
                        normal_gr=normal_gr,
                        normal_raw=normal_raw,
                        normal_cov=normal_cov)

TS <- sort(TS)

saveWithID(TS, "sample_specific_variants_granges", SAMid,
           save_dir = file.path(save_dir),
           compress = FALSE)

VariantTools::cgpGr2vcf(GR=TS, sample_id=SAMid,
          filename=file.path(save_dir,
            paste(SAMid,"_sample_specific_variants.vcf", sep = "")),
          project=project)

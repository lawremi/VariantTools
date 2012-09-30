library(VariantTools)
##library(RNASeqGenie)


cmd_args = commandArgs(trailingOnly=TRUE);
bam_file <- cmd_args[1];
read_length=cmd_args[2];
genome<- cmd_args[3];
genome_dir <- cmd_args[4];
save_dir <- cmd_args[5];
SAMid <- cmd_args[6];
project<-cmd_args[7];

if(file.exists(save_dir)){
  stop("Error: Save directory already exists")
}else{
  dir.create(save_dir)
  dir.create(file.path(save_dir, "results"))
}
           
qual <- detectQualityType(bam_file)
if(qual=="illumina"){
  var<- callVariantsP(bam=bam_file, genome=genome, genome_dir=genome_dir,mapq=13, bqual=56, force=TRUE,with_qual=T, read_length=read_length)
}else if(qual=="sanger"){
  var<- callVariantsP(bam=bam_file, genome=genome, genome_dir=genome_dir,mapq=13, bqual=23,with_qual=T, read_length=read_length)
}else{
  stop("Error: Non-standard quality in BAM")
}
saveWithID(var$raw_granges, "variants_granges", SAMid,
           save_dir = file.path(save_dir, "results"),
           compress = FALSE)
saveWithID(var$filtered_granges, "filtered_variants_granges", SAMid,
           save_dir = file.path(save_dir,"results"),
           compress = FALSE)

VariantTools::cgpGr2vcf(GR=var$filtered_granges, sample_id=SAMid,
          filename=file.path(save_dir, "results",
            paste(SAMid,".vcf", sep = "")),
          project=project)

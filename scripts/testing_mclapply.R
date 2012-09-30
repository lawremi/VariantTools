

genome_gr <- getChrGRangesFromGmap(genome = genome, genome_dir = genome_dir)
chr_ids <- as.character(seqnames(genome_gr))
callVariantsP<-function(merged_bam,
                        read_length=-1,
                        genome=NULL,
                        genome_dir=NULL,
                        with_qual=FALSE,
                        bqual = 0,
                        mapq=0,
                        mc.cores = detectCores()){
  
  gr_list <- mclapply(chr_ids, mc.cores=mc.cores, function(x){
    chr_tally_gr <-tally2GR(bamfiles= merged_bam,
                            genome =genome,
                            genome_dir=genome_dir,
                            regions = genome_gr,
                            breaks=NULL,
                            bqual_thresh=56,
                            map_qual=13,
                            mc.cores = 1,
                            chr_ids=x)
    
    seqlevels(chr_tally_gr) =seqlevels(genome_gr)
    seqlengths(chr_tally_gr) =seqlengths(genome_gr)
    
    chr_variants_filtered <- variantFilter(granges=chr_tally_gr, useQual=TRUE)
  })
  tally_gr <- Reduce(f=merge_var, gr_list)
  return(tally_gr)
}
merge_var <- function(a, b){
  res <- list()
  res$raw_granges<- c(a$raw_granges, b$raw_granges)
  res$filtered_granges <- c(a$filtered_granges, b$filtered_granges)
  res$N_nucleotide_rejects <- sum( a$N_nucleotide_rejects, b$N_nucleotide_rejects, na.rm=T)
  res$low_count_rejects<-sum(a$low_count_rejects, b$low_count_rejects, na.rm=T)
  res$frequency_rejects<- sum(a$frequency_rejects, b$frequency_rejects, na.rm=T)
  res$fisher_test_rejects_2by2<-sum(a$fisher_test_rejects_2by2, b$fisher_test_rejects_2by2,na.rm=T)
  res$fisher_test_rejects_3by2<-sum(a$fisher_test_rejects_3by2, b$fisher_test_rejects_3by2, na.rm=T)
  return(res)
}

merged_bam <- "/gne/home/degenhj2/svn/Packages/VariantTools/tests/merged/SRC111111.TESTPREPENDSTR.merged.uniques.bam"
genome = "hg19"
genome_dir <- "/gnet/is2/data/bioinfo/gmap/data/genomes/"
genome_gr <- getChrGRangesFromGmap(genome = genome, genome_dir = genome_dir)
chr_ids <- as.character(seqnames(genome_gr))

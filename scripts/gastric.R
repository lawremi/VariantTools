path <- "/gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS169/"
sams <- dir(path)
genome_dir = "/gne/research/workspace/degenhj2/"
genome = "hg19_IGIS21"
post = ".analyzed.bam"
sams <- as.list(sams)
calls <- mclapply(sams, mc.cores = detectCores(), function(sam){
  bam_file <- file.path(path,sam,"bams/", paste(sam, post, sep = ""))
  gr <- tally2GR(bamfiles=bam_file, genome=genome,genome_dir=genome_dir,chr_gr=gene, variant_strand=1, map_qual=13, bqual_thresh=23)
  gr <-variantFilter(gr, useQual=TRUE)
  gr <- gr$filtered_granges
  return(gr)
}
)

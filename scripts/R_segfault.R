library(VariantTools)
library(parallel)
bam <- "/gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS121/SAM587376/bams/SAM587376.analyzed.bam"
genome = "hg19_IGIS21"
genome_dir = "/gne/research/workspace/degenhj2/"
bam = file.path(system.file(package = "VariantTools", mustWork=TRUE),"extdata/merged/SRC111111.TESTPREPENDSTR.merged.uniques.bam")
genome = "hg19_ucsc"
genome_dir = "/gnet/is2/research/data/bioinfo/gmap/data/genomes"
bqual = 23
mapq=13
mc.cores = detectCores()
force_qual=TRUE
read_length=75
bin_fraction=0.2
genome <- GmapGenome(genome=genome, directory=genome_dir)
genome_gr <- as(seqinfo(genome), "GenomicRanges")

chr_inds <-seq_len(length(genome_gr))
qual_type <-detectQualityType(bam, 500)
breaks <- c(0, round(bin_fraction*read_length), round((1-bin_fraction)*read_length), read_length)
read_pos <- paste('cycleCount.', breaks[2], ".", breaks[3], sep ="")
useQual =TRUE
gr_list <- mclapply(chr_inds, mc.cores=mc.cores, function(x){
  chr_info <- genome_gr[x]
  bf <- GmapBamReader(bam)
  chr_tally_gr <-bam_tally(bf,
                           genome,
                           which=chr_info,
                           cycle_breaks= breaks,
                           minimum_mapq=mapq,
                           variant_strand=1,
                           high_quality_cutoff=bqual)
  chr_variants_filtered <- variantFilter(granges=chr_tally_gr, useQual=useQual, read_pos=read_pos)
  
})
var <- VariantTools:::merge_var(gr_list)

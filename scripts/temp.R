gsnapChromosomes <- function(genome, genome_dir) {
  read.table(pipe(paste("get-genome -L -d", genome, "-D", genome_dir)))[,1]
}


tally2GR<- function(bamfiles,
                    genome = "hg19_ucsc",
                    genome_dir = "/gnet/is2/data/bioinfo/gmap/data/genomes",
                    chr_ids = gsnapChromosomes(genome, genome_dir),
                    regions,
                    variant_strand=1,
                    min_cov =1,
                    breaks=NULL,
                    bqual_thresh=56,
                    map_qual=13)
{
  chr_ids <- as.list(as.character(chr_ids))
  has_regions <- !missing(regions)
  list_of_gr <- mclapply(chr_ids, function(chr_name){
    tally <- pipe(paste("bam_tally -B 0 -C -Q -q", map_qual, " -X ", variant_strand, " -n ",
                        min_cov, " -T -d", genome, "-D", genome_dir, bamfiles,
                        paste("'", chr_name, ":' 2> /dev/null", sep = "")))
    tab <- read.table(tally, colClasses = c("character", "integer", "integer",
                               "character"), sep = "\t",
                      col.names = c("chrom", "position", "count", "cycles"))
    if(dim(tab)[1] <1){
      gr <- GRanges()
    } else {      
      counts <- strsplit(tab[[4]], " ", fixed=TRUE)
      
      counts_part <- PartitioningByWidth(elementLengths(counts))
      
      chr <- rep(tab[[1]], width(counts_part))
      pos <- rep(tab[[2]], width(counts_part))
      count.total <- rep(tab[[3]], width(counts_part))
      ## we only get duplicated locations if there is overlap across strands
      ## thus, we force the duplicated locations onto different strands
      
      location <- paste(tab[[1]], tab[[2]], sep = ":")
      location_dup <- duplicated(location)
      strand <- rep("*", length(location))
      
      strand[match(location[location_dup], location)] <- "+"
      strand[location_dup] <- "-"
      strand <- rep(strand, width(counts_part))
      location <- rep(location, width(counts_part))
      counts_flat <- unlist(counts, use.names=FALSE)
      bases <- sub("(\\D).*", "\\1", counts_flat)
      ref_rows <- start(counts_part)
      ref <- rep(bases[ref_rows], width(counts_part))
      zero_count <- !grepl("(", counts_flat, fixed=TRUE)
### FIXME: obscene hack to simplify code
      counts_flat[zero_count] <- paste(counts_flat[zero_count], "(0@NA)", sep="")
      cycles <- strsplit(sub(".*\\((.*?)\\).*", "\\1", counts_flat), ",", 
                         fixed=TRUE)
      ncycles <- elementLengths(cycles)
      cycles_flat <- unlist(cycles, use.names=FALSE)
      cycles_counts <- unlist(strsplit(cycles_flat, "@", fixed=TRUE), 
                              use.names=FALSE)
      cycles_mat <- matrix(suppressWarnings(as.integer(cycles_counts)), nrow=2)
      ## insane code for counting up reads by strand
      neg_cycle <- cycles_mat[2,] < 0
      coord_ind <- rep(seq_along(cycles), ncycles)
      strand_key <- paste(coord_ind, neg_cycle)
      strand_key_uniq <- !duplicated(strand_key)
      strand_key <- factor(strand_key, strand_key[strand_key_uniq])
      neg_cycle_uniq <- neg_cycle[strand_key_uniq]
      coord_ind_uniq <- coord_ind[strand_key_uniq]
      strand_count <- rowsum(cycles_mat[1,], strand_key)
      count.pos <- count.neg <- rep(0L, length(cycles))
      ## have to use which() here to drop the NA values
      count.pos[coord_ind_uniq[which(!neg_cycle_uniq)]] <-
        strand_count[which(!neg_cycle_uniq)]
      count.neg[coord_ind_uniq[which(neg_cycle_uniq)]] <-
        strand_count[which(neg_cycle_uniq)]
      message("working on breaks") 
      if (!is.null(breaks)) {
### NOTE: overuse of rep() here might lead to over-long vectors            
        cycle_bins <- cut(rep(abs(cycles_mat[2,]), cycles_mat[1,]), breaks)      
        cycle_i <- factor(rep(rep(seq(length(cycles)), ncycles), cycles_mat[1,]),
                          seq(length(cycles)))
        cycle_tab <- table(cycle_i, cycle_bins)
        colnames(cycle_tab) <- paste(head(breaks, -1), tail(breaks, -1),
                                     sep = ".")
        ref_cycle_tab <- cycle_tab[rep(ref_rows, width(counts_part)),]
        colnames(ref_cycle_tab) <- paste(colnames(cycle_tab), "ref", sep = ".")
        cycle_tab <- cbind(cycle_tab, ref_cycle_tab)
      } else  cycle_tab <- matrix(nrow = length(pos), ncol = 0)
      message("completed breaks")
      base_counts <- as.integer(sub("\\D(\\d+).*", "\\1", counts_flat))
      ncycles[zero_count] <- 0L # otherwise 1 due to our obscene hack above
      ncycles.ref <- rep(ncycles[ref_rows], width(counts_part))
      count.ref <- rep(base_counts[ref_rows], width(counts_part))
      count.pos[zero_count] <- 0L
      count.pos.ref <- rep(count.pos[ref_rows], width(counts_part))
      count.neg[zero_count] <- 0L
      count.neg.ref <- rep(count.neg[ref_rows], width(counts_part))
      quals <- strsplit(sub(".*\\)\\((.*?)\\).*", "\\1", counts_flat), ",", 
                        fixed=TRUE)
      nquals <- elementLengths(quals)

      count_above_thresh <- lapply(quals, function(x){
        mat <-matrix(unlist(strsplit(x, "Q", fixed=TRUE), 
                      use.names=FALSE), nrow=2)
        logical <- (as.numeric(mat[2,]))>bqual_thresh
        sum(as.numeric(mat[1,logical]), na.rm=T)
      })
      high_qual <- unlist(count_above_thresh)      
      gr <- GRanges(chr, IRanges(pos, width=1L), strand, location, 
                    ref = DNAStringSet(ref), read = DNAStringSet(bases),
                    ncycles, ncycles.ref,
                    count = base_counts, count.ref, count.total, high.quality = high_qual,
                    count.pos, count.pos.ref, count.neg, count.neg.ref,
                    cycleCount = unclass(cycle_tab))
      
      ##adding ref counts back in as counts per break segment
      ##there seems to be a bug somewhere above that is causing the ref positions to be listed twice some times. Removing these for now
      gr <- gr[as.character(values(gr)$ref) != as.character(values(gr)$read)]

      if (has_regions) {
        region_ol <- findOverlaps(gr, regions)
        region_strand <- as.vector(strand(regions))[subjectHits(region_ol)]
        strand(gr) <- region_strand
        rc <- region_strand == "-"
        values(gr)$ref[rc] <- reverseComplement(values(gr)$ref[rc])
        values(gr)$read[rc] <- reverseComplement(values(gr)$read[rc])
      }
      values(gr)$location <- paste(values(gr)$location, strand(gr), sep = ":")
    }
    message(paste("finished chr ", chr_name))
    seqlevels(gr) <- unlist(chr_ids)
    gr
  })
  GR_full <- do.call(c, list_of_gr)
  return(GR_full)
}




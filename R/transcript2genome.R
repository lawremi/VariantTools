trans2genome <-function(tx_name, cds_db, pos, type=c("codon", "nuc")){
  strand=suppressWarnings(as.character(unlist(strand(cds_db[tx_name])))[1])
  if(type=="codon"){
    pos=.codon2trans(pos)
  }
  genome_pos <-transcriptLocs2refLocs(list(as.integer(pos)),
                         start(cds_db[tx_name]),
                         end(cds_db[tx_name]),
                         strand=strand
                                      )
  if(type=="codon"){
    if(strand =="-"){
      gr <- GRanges(seqnames=suppressWarnings(as.character(unlist(seqnames(cds_db[tx_name])))[1]),
                    IRanges(start=(unlist(genome_pos)-2), width=3))
    }
    if(strand=="+"){
      gr <- GRanges(seqnames=suppressWarnings(as.character(unlist(seqnames(cds_db[tx_name])))[1]),
                    IRanges(start=(unlist(genome_pos)), width=3))
    }
  }else{
    gr <- GRanges(seqnames=suppressWarnings(as.character(unlist(seqnames(cds_db[tx_name])))[1]),
                  IRanges(start=(unlist(genome_pos)), width=1))
  }
  
  return(gr)
}

.codon2trans<-function(codon){
  pos <- ((codon-1)*3)+1
  return(pos)
}

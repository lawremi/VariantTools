rnaValRate <- function(gr, rna, cov, type = "rate",min=0){
  values(gr)$norm_ref <- as.numeric(NA)
  values(gr)$norm_alt <- as.numeric(NA)
  over2 <- which(paste(values(gr)$location,values(gr)$read) %in% paste(values(rna)$location, values(rna)$read))
  over <- which(paste(values(rna)$location,values(rna)$read) %in% paste(values(gr)$location, values(gr)$read))
  values(gr[over2])$norm_ref <- values(rna[over])$count.ref
  values(gr[over2])$norm_alt <- values(rna[over])$count
  values(gr[-over2])$norm_alt <- 0
  gr<-gr[values(gr)$count.total>min]
  grl <- split(gr, as.factor(seqnames(gr)))
  names <- as.list(names(grl))
  n <- lapply(names, function(x) { as.numeric(cov[[x]][ranges(grl[[x]])]) })
  n <- unlist(n)
  gr <- unlist(grl)
  values(gr)$norm_cov<-n
  gr[values(gr)$norm_cov >20]->high
  if(length(high)>=0){
    val_rate <- 1- length(high[values(high)$norm_alt==0])/length(high)
    cat("RNA validation rate: ", val_rate, "\n", sep = "")
    if(type == "rate"){
      return(c(val_rate, length(high), length(high[values(high)$norm_alt==0])))
    }else{
      return(gr)
    }
  }else{
    return(NULL)
  }
}

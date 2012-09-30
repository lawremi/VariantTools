detectReadLength <- function(bamfile, nrows=10){
  bf <- BamFile(file = bamfile, yieldSize=nrows)
  param <- ScanBamParam(what="seq")
  bam <- scanBam(bf, param=param)
  reads <- as.character(bam[[1]]$seq)
  length <- length(reads)
  if(length ==0){
    stop("no reads in bam file")
  }
  read_len <- mean(nchar(reads))
  return(read_len)
}
  

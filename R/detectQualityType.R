###############################
## Detect the type of quality
## scores encoded in a bam file
###############################

detectQualityType <- function(bamfile, nrows=5000){
  qual_type <- array()
  bf <- BamFile(file = bamfile, yieldSize=nrows)
  param <- ScanBamParam(what="qual")
  bam <- scanBam(bf, param=param)
  quals <- as.character(bam[[1]]$qual)
  length <- length(quals)
  if(length ==0){
    stop("no reads in bam file")
  }
  max_vec <- array()
  min_vec <- array()
  for(i in seq_along(length)){
    max_vec[i]<-max(as.integer(charToRaw(quals[i])))
    min_vec[i]<-min(as.integer(charToRaw(quals[i])))
  }
  min <- min(min_vec)
  max <- max(max_vec)
  if(max >74 & min >=59){
    qual_type[1]="illumina"
  }else if(max <=74 & min >= 33){
    qual_type[1]="sanger"
  }else if(min<33 & max <= 40){
    qual_type[1]="phred"
  }else{
    qual_type[1]="fubar"
  }
  return(qual_type)
}

checkQualityType <- function(bam, expected) {
  detected <- detectQualityType(bam, 500)
  if(detected != expected)
    {
      stop("Non-standard qualities in BAM. Quality type returned: ",
           detected, ", expected: ", expected)
    }
  invisible(TRUE)
}

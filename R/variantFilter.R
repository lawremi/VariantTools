##  Function to filter nextgen seq data to detect variants
##

## ==========================================================
## Original notes from Jeremiah
## ==========================================================

## The function first checks and either performs Fisher's Exact
## Test(s) and appends the pvalues onto the granges, and saves the
## granges, or it saves the raw granges and removes strands with 'N'
## nucleotide variant or reference occurences and strands with lower
## than 2 counts per variant, depending on the input parameter. It
## next filters out all granges with Fisher's Test pvals for either
## test less than the designated parameter. It then filters out all
## GRange data entries that have lower coverage as well as lower
## variant occurence frequencies than those designated by the
## parameters.

## ==========================================================
## What we think the function actually does
## ==========================================================
##
## 1) Considers only positions where ALT != REF
## 2) Discards the variants with N in ALT or REF
## 3) Discards the variants with low ALT read count (< 'readCount=2')
## 4) Discards the variants with low ALT cycle count (< 'cycleCount=2')
## 5) If 'lrt' but not 'lr', then perform binomial test, with 'error' prob
##    and upper-tail p-value cutoff ('lrt_p')
## 6) If 'lr' (preferred?), then calculate P(D|p=0.2) / P(D|p=error)
##    and check if it is higher than one. If 'useQual', then this uses
##    the high quality fraction, instead of total.
## 7) If 'fet_2x2' it checks for strand bias using FET, < pval = 0.05
## 8) If 'read_pos' is given as a column name, presumably matching a
##    cycleCount bin column, it checks that the column is > 0. Typically,
##    this will be a middle bin (since the ends are not as trusted).

variantFilter <- function(granges,
                          useQual=FALSE,
                          pval = 0.05,
                          readCount=2,
                          cycleCount=1,
                          lrt = FALSE,
                          lrt_p = 0.01,
                          lr=TRUE,
                          lower = 0.2,
                          error=1/1000,
                          fet_2x2 = TRUE,
                          read_pos=NULL
                          ) {
############################################################################
  
  output_list_names <- c("raw_granges", "filtered_granges", "N_nucleotide_rejects",
                         "low_count_rejects", "frequency_rejects", "fisher_test_rejects_2by2",
                         "read_pos_rejects")
  output <- list(granges, GRanges(), 0, 0,0,0,0)
  names(output) <- output_list_names
  
  if(length(granges) == 0){
    ##message("no variants found!")
    output <- list(granges, granges, 0, 0,0,0,0)
    names(output) <- output_list_names
    return(output)
  }
  granges <- granges[as.character(values(granges)$ref) != as.character(values(granges)$alt)]
  
############################################################################
  ##removing N nucleotides and low counts of alternate reads
  Ngranges <- granges[as.character(values(granges)$alt) !='N']
  Ngranges <- Ngranges[as.character(values(Ngranges)$ref) !='N'] 
  Nrejects <- length(granges) - length(Ngranges)
  lowcountgranges <- Ngranges[values(Ngranges)$count >= readCount]
  if(!is.null(cycleCount)){
    lowcountgranges <- lowcountgranges[values(lowcountgranges)$ncycles > cycleCount]
  }
  lowcountrejects <- length(Ngranges) - length(lowcountgranges)
  ##if no more positions return
  if(length(lowcountgranges)==0){
    output[[2]]<-lowcountgranges
    output[[3]] <- Nrejects
    output[[4]]<-lowcountrejects
    return(output)
  }
  
############################################################################
  if(lrt ==TRUE){
    ##message("LRT on frequency...")
    res <-pbinom((values(lowcountgranges)$count-1),
                 values(lowcountgranges)$count.total,
                 prob=error,
                 lower.tail = FALSE
                 )
    freq_granges <- lowcountgranges[res < lrt_p]
  }
  if(lr ==TRUE){
    freq <- freq(error, lower, 1)
    if(useQual==TRUE){
      lowcountgranges <- lowcountgranges[(values(lowcountgranges)$high.quality+values(lowcountgranges)$high.quality.ref) !=0]
      freq_granges <- lowcountgranges[(values(lowcountgranges)$high.quality /
                                       (values(lowcountgranges)$high.quality+values(lowcountgranges)$high.quality.ref)) >= freq]
    }else{
     freq_granges <- lowcountgranges[(values(lowcountgranges)$count/values(lowcountgranges)$count.total) >= freq]
   }
  }
  freqrejects <- length(lowcountgranges) - length(freq_granges)
  if(length(freq_granges)==0){
    output[[2]]<-freq_granges;
    output[[3]] <- Nrejects
    output[[4]]<-lowcountrejects
    output[[5]]<- freqrejects
    return(output)
  }
############################################################################
  ##Performing the Fisher's Test Filters
  if(fet_2x2 == TRUE){
    ## message("Filtering granges by Fishtest pvals...")
    two_by_2_pvals <- list()
    x1 <- fisher_p(values(freq_granges)$count.pos.ref,
                   values(freq_granges)$count.neg.ref,
                   (values(freq_granges)$count.pos.ref+values(freq_granges)$count.pos),
                   (values(freq_granges)$count.neg.ref+values(freq_granges)$count.neg))
    two_by_two_filter_gr <- freq_granges[x1 >= pval]
    fishers_2by2_rejects <- length(freq_granges) - length(two_by_two_filter_gr)
    if(length(two_by_two_filter_gr)==0){
      output[[2]]<-two_by_two_filter_gr
      output[[3]] <- Nrejects
      output[[6]]<-lowcountrejects
      output[[5]]<-freqrejects
      output[[6]]<-fishers_2by2_rejects
      return(output)
    }
  }else{
    two_by_two_filter_gr <- freq_granges
  }
  if(!is.null(read_pos)){
    output_granges<-two_by_two_filter_gr[values(two_by_two_filter_gr)[[read_pos]] >0]
    read_pos_rejects <- length(two_by_two_filter_gr) - length(output_granges)
  }else{
    output_granges<-two_by_two_filter_gr
  }
############################################################################

  output[['filtered_granges']] <- output_granges 
  output[[3]] <- Nrejects
  output[[4]] <- lowcountrejects
  output[[5]] <- freqrejects
  if (fet_2x2 == TRUE) {
    output[[6]] <- fishers_2by2_rejects
  } else {
    output[[6]] <- 0
  }
  if (!is.null(read_pos)) {
    output[[7]] <- read_pos_rejects
  } else {
    output[[7]] <- 0
  }
  ## message("done!")
  return(output)
}

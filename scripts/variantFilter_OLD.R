
##  Function to filter nextgen seq data to detect variants
##
## The function first checks and either performs Fisher's Exact Test(s) and appends the pvalues onto the granges,
## and saves the granges, or it saves the raw granges and removes  strands with 'N' nucleotide variant or reference
## occurences and strands with lower than 2 counts per variant, depending on the input parameter. It next filters
## out all granges with Fisher's Test pvals for either test less than the designated parameter. It then filters out all
## GRange data entries that have lower coverage as well as lower variant occurence frequencies than those
## designated by the parameters. Finally, if the parameter dictates that the function should not append pvals,
## it finally performs 2x2 and 3x2 Fisher's Exact Tests tofilter for independence between ref/variant strands and
## between positive and negative strands. 

variantFilter <- function(granges,
                          useQual=FALSE,
                          pval = 0.05,
                          lrt = FALSE,
                          lrt.p = 0.01,
                          lr=TRUE,
                          lower = 0.2,
                          error=1/1000,
                          fet_3x2 = FALSE,
                          fet_2x2 = TRUE,
                          append = FALSE) {
############################################################################

  output_list_names <- c("raw_granges", "filtered_granges", "N_nucleotide_rejects",
                         "low_count_rejects", "frequency_rejects", "fisher_test_rejects_2by2",
                         "fisher_test_rejects_3by2")
  output <- list(GRanges(), GRanges(), 0, 0,0,0,0)
  names(output) <- output_list_names
  
  if(length(granges) == 0){
    ##message("no variants found!")
    output <- list(granges, granges, 0, 0,0,0,0)
    names(output) <- output_list_names
    return(output)
  }
  granges <- granges[as.character(values(granges)$ref) != as.character(values(granges)$read)]
  if (append == TRUE) {
    ##saving Fisher's Test pvalues onto the raw granges
    if(fet_2x2 == TRUE){
   ##   message("Saving Fisher's Tests pvals to granges...") 
      two_by_2_pvals <- list()
      temp_mat1 <- cbind(values(granges)$count.pos.ref,
                         values(granges)$count.pos,
                         values(granges)$count.neg.ref,
                         values(granges)$count.neg)
      x1 <- vector(mode = 'numeric', length = nrow(temp_mat1))
      for(j in seq_len(nrow(temp_mat1))) {
        if (sum(as.numeric(temp_mat1[j,])) < 90000) {
          x1[j] <- fisher.test(matrix(temp_mat1[j,], ncol = 2, byrow = FALSE), workspace = 1000000)$p.val
        } else {
          x1[j] <- fisher.test(matrix(temp_mat1[j,], ncol = 2, byrow = FALSE), workspace = 1000000, simulate.p.value = TRUE)$p.val
        }
      } 
      two_by_2_pvals <- as.list(x1)
      elementMetadata(granges)$fishers2by2.pval <- x1
    } else if (fet_3x2 == TRUE) {
      warning("Hard-coded to expect tally2GR was given three bins. Generalize this")
      three_by_2_pvals <- list()
      temp_mat <- cbind(values(granges)[,16],
                        values(granges)[,13],
                        values(granges)[,17],
                        values(granges)[,14],
                        values(granges)[,18],
                        values(granges)[,15])
      x <- vector(mode = 'numeric', length = nrow(temp_mat))
      for(j in seq_len(nrow(temp_mat))) {
        if (sum(as.numeric(temp_mat[j,])) < 90000) {
          x[j] <- fisher.test(matrix(temp_mat[j,], ncol = 3, byrow = FALSE), workspace = 1000000)$p.val
        } else {
          x[j] <- fisher.test(matrix(temp_mat[j,], ncol = 3, byrow = FALSE), workspace = 1000000, simulate.p.value = TRUE)$p.val
        }
      }
      three_by_2_pvals <- as.list(x)
      elementMetadata(granges)$fishers3by2.pval <- x
    }
    granges1 <- granges
    output[[1]] <- granges1
    if (length(granges1) == 0)
      return(output)
    
    ##Performing the Fisher's Test Filters
    ##message("Filtering granges by Fishtest pvals...")
    if(fet_2x2 == TRUE){
      two_by_two_filter_gr <- granges[values(granges)$fishers2by2.pval >= pval]
      fishers_2by2_rejects <- length(granges) - length(two_by_two_filter_gr)
    }else{
      two_by_two_filter_gr <- granges
    }
    if(fet_3x2 == TRUE){
      three_by_two_filter_gr <- two_by_two_filter_gr[values(two_by_two_filter_gr)$fishers3by2.pval >= pval]
      fishers_3by2_rejects <- length(two_by_two_filter_gr) - length(three_by_two_filter_gr)
      granges <- three_by_two_filter_gr
    }
  }
   output[['raw_granges']] <- granges
############################################################################
  ##removing N nucleotides and low counts of alternate reads
  Ngranges <- granges[as.character(values(granges)$read) !='N']
  Ngranges <- Ngranges[as.character(values(Ngranges)$ref) !='N'] 
  Nrejects <- length(granges) - length(Ngranges)
  lowcountgranges <- Ngranges[values(Ngranges)$count >= 2]
  lowcountrejects <- length(Ngranges) - length(lowcountgranges)
  ##if no more positions return
  if(lowcountrejects==length(Ngranges)){
    output[[2]]<-lowcountgranges;
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
    freq_granges <- lowcountgranges[res < lrt.p]
  }
  if(lr ==TRUE){
    freq <- freq(error, lower, 1)
    if(useQual==TRUE){
      lowcountgranges <- lowcountgranges[(values(lowcountgranges)$high.quality+values(lowcountgranges)$high.quality.ref) !=0]
      freq_granges <- lowcountgranges[(values(lowcountgranges)$high.quality/(values(lowcountgranges)$high.quality+values(lowcountgranges)$high.quality.ref)) >= freq]
    }else{
     freq_granges <- lowcountgranges[(values(lowcountgranges)$count/values(lowcountgranges)$count.total) >= freq]
   }
  }
     freqrejects <- length(lowcountgranges) - length(freq_granges)
     if(freqrejects==length(lowcountgranges)){
       output[[2]]<-freq_granges;
       output[[3]] <- Nrejects
       output[[4]]<-lowcountrejects
         output[[5]]<- freqrejects
       return(output)
     }
        if(append == TRUE) {
    output_granges <- freq_granges
  }
############################################################################
  if (!append) {
    ##Performing the Fisher's Test Filters
    if(fet_2x2 == TRUE){
     ## message("Filtering granges by Fishtest pvals...")
      two_by_2_pvals <- list()
      temp_mat1 <- cbind(values(freq_granges)$count.pos.ref,
                         values(freq_granges)$count.pos,
                         values(freq_granges)$count.neg.ref,
                         values(freq_granges)$count.neg)
      x1 <- vector(mode = 'numeric', length = nrow(temp_mat1))
      for(j in seq_len(nrow(temp_mat1))) {
        if (sum(as.numeric(temp_mat1[j,])) < 90000) {
          x1[j] <- fisher.test(matrix(temp_mat1[j,], ncol = 2, byrow = F), workspace = 1000000)$p.val
        } else {
          x1[j] <- fisher.test(matrix(temp_mat1[j,], ncol = 2, byrow = F), workspace = 1000000, simulate.p.value = TRUE)$p.val
        }
      }
      index <- which(x1 >= pval)
      two_by_two_filter_gr <- freq_granges[index]
      fishers_2by2_rejects <- length(freq_granges) - length(two_by_two_filter_gr)
      if(fishers_2by2_rejects==length(freq_granges)){
        output[[2]]<-lowcountgranges;
        output[[3]] <- Nrejects
        output[[6]]<-lowcountrejects
          output[[5]]<-freqrejects
        output[[6]]<-fishers_2by2_rejects
        return(output)
      }
    }else{
        two_by_two_filter_gr <- freq_granges
    }
    if(fet_3x2 == TRUE){
      warning("Hard-coded to expect tally2GR was given three bins. Generalize this")
      three_by_2_pvals <- list()
      temp_mat <- cbind(values(two_by_two_filter_gr)[,16],
                        values(two_by_two_filter_gr)[,13],
                        values(two_by_two_filter_gr)[,17],
                        values(two_by_two_filter_gr)[,14],
                        values(two_by_two_filter_gr)[,18],
                        values(two_by_two_filter_gr)[,15])
      x <- vector(mode = 'numeric', length = nrow(temp_mat))
      for(i in seq_len(nrow(temp_mat))) {
        if (sum(as.numeric(temp_mat[i,])) < 90000) {
          x[i] <- fisher.test(matrix(temp_mat[i,], ncol = 3, byrow = FALSE), workspace = 1000000)$p.val
        } else {
          x[i] <- fisher.test(matrix(temp_mat[i,], ncol = 3, byrow = FALSE), workspace = 1000000, simulate.p.value = TRUE)$p.val
        }
      }
      index2 <- which(x >= pval)
      output_granges <- two_by_two_filter_gr[index2]
      fishers_3by2_rejects <- length(two_by_two_filter_gr) - length(output_granges)
    }else{
      output_granges <- two_by_two_filter_gr
    }
  }
############################################################################

  output[['filtered_granges']] <- output_granges
  if (length(output_granges) == 0)
    return(output)    

  ##message("outputing result...")
  output[[3]] <- Nrejects
  output[[4]] <- lowcountrejects
  output[[5]] <- freqrejects
  if (fet_2x2 == TRUE) {
    output[[6]] <- fishers_2by2_rejects
  } else {
    output[[6]] <- 0
  }
  if (fet_3x2 == TRUE) {
    output[[7]] <- fishers_3by2_rejects
  } else {
    output[[7]] <- 0
  }
 ## message("done!")
  return(output)
}


freq <- function(p0,p1, n){
  num <- (log(1-p0)-log(1-p1))
  denom <- (log(p0)-log(p1)+log(1-p1)-log(1-p0))
  rat <- num/denom
  ret <- -n*rat
  return(ret)
}

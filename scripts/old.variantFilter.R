##'  Function to filter raw varaint positions found in a bam file
##'
##' The function first checks and either performs Fisher's Exact Test(s) and appends the pvalues onto the granges,
##' and saves the granges, or it saves the raw granges and removes  strands with 'N' nucleotide variant or reference
##' occurences and strands with lower than 2 counts per variant, depending on the input parameter. It next filters
##' out all granges with Fisher's Test pvals for either test less than the designated parameter. It then filters out all
##' GRange data entries that have lower coverage as well as lower variant occurence frequencies than those
##' designated by the parameters. Finally, if the parameter dictates that the function should not append pvals,
##' it finally performs 2x2 and 3x2 Fisher's Exact Tests tofilter for independence between ref/variant strands and
##' between positive and negative strands. 
##' 
##' @title variantFilter
##' @param granges Grange object of variance appropriately formated from tally.R of the ExomeSeq package
##' @param cov_val A coverage threshold
##' @param cov An rle vector denoting the coverage from the bam file from the sequencing lane
##' @param pval The pval threshold for the Fisher's Exact Tests
##' @param freq The frequency threshold for variance
##' @param prepend_str An expression or ID to prepend to the grange file title
##' @param save_dir A directory into which the files will be saved
##' @param src_id The lane's SRC identification
##' @param fet_3x2 A logical parameter dictating whether the function should perform a 3x2 Fisher's Exact Test
##' on the strands.
##' @param fet_2x2 A logical parameter dictating whether the function should perform a 2x2 Fisher's Exact Test
##' on the strands.
##' @param append A logical parameter that dictates whether the function appends the Fisher's Exact Test
##' pvalues onto the granges or not.
##' @return A list containing the raw Granges, filtered GRange data set, the number of N nucleotide rejects, the
##' rejects with fewer than 2 counts, the number of coverage rejects, the number of frequency rejects, the
##' number of rejects from the 2x2 Fisher's, and the number of rejects from the 3x2 Fisher's in that order. Also
##' the vcf data and the raw granges are printed in separate files.
##' @author Jeremiah Degenhardt and Jason Young 
##' @export
variantFilter <- function(granges,
                          pval = 0.05,
                          binom.p = 0.01,
                          fet_3x2 = FALSE,
                          fet_2x2 = TRUE,
                          append = FALSE) {
############################################################################

  if (append  == TRUE) {
    ##saving Fisher's Test pvalues onto the raw granges
    if(fet_2x2 == TRUE){
      message("Saving Fisher's Tests pvals to granges...") 
      two_by_2_pvals <- list()
      temp_mat1 <- cbind(values(granges)$count.pos.ref,
                         values(granges)$count.pos,
                         values(granges)$count.neg.ref,
                         values(granges)$count.neg)
      x1 <- vector(mode = 'numeric', length = dim(temp_mat1)[1])
      for(j in 1:dim(temp_mat1)[1]) {
        if (sum(as.numeric(temp_mat1[j,])) < 90000) {
          x1[j] <- fisher.test(matrix(temp_mat1[j,], ncol = 2, byrow = F), workspace = 1000000)$p.val
        } else {
          x1[j] <- fisher.test(matrix(temp_mat1[j,], ncol = 2, byrow = F), workspace = 1000000, simulate.p.value = TRUE)$p.val
        }
      } 
      two_by_2_pvals <- as.list(x1)
      elementMetadata(granges)$fishers2by2.pval <- x1
    }else{
    }
    if(fet_3x2 == TRUE){
      warning("Hard-coded to expect tally2GR was given three bins. Generalize this")
      three_by_2_pvals <- list()
      temp_mat <- cbind(
                        values(granges)[,16],
                        values(granges)[,13],
                        values(granges)[,17],
                        values(granges)[,14],
                        values(granges)[,18],
                        values(granges)[,15]
                        )
      x <- vector(mode = 'numeric', length = dim(temp_mat)[1])
      for(j in 1:dim(temp_mat)[1]) {
        if (sum(as.numeric(temp_mat[j,])) < 90000) {
          x[j] <- fisher.test(matrix(temp_mat[j,], ncol = 3, byrow = F), workspace = 1000000)$p.val
        } else {
          x[j] <- fisher.test(matrix(temp_mat[j,], ncol = 3, byrow = F), workspace = 1000000, simulate.p.value = TRUE)$p.val
        }
      }
      three_by_2_pvals <- as.list(x)
      elementMetadata(granges)$fishers3by2.pval <- x
      ##loginfo("...done")
    }else{
    }
    granges1 <- granges
    ##Performing the Fisher's Test Filters
    message("Filtering granges by Fishtest pvals...")
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
    ##loginfo("...done") 
  } else {
  
}
############################################################################
  ##removing N nucleotides and low counts of alternate reads
  message("filtering granges for N nucleotides and variant counts less than 2 per variant...")
  Ngranges <- granges[as.character(values(granges)$read) !='N']
  Ngranges <- Ngranges[as.character(values(Ngranges)$ref) !='N']
  Nrejects <- length(granges) - length(Ngranges)
  lowcountgranges <- Ngranges[values(Ngranges)$count >= 2]
  lowcountrejects <- length(Ngranges) - length(lowcountgranges)
  ##loginfo("...done")
############################################################################

### update this version
  message("LRT on frequency...")
   res <-   pbinom(
                   (values(lowcountgranges)$count-1),
                   values(lowcountgranges)$count.total,
                   prob=1/1000,
                   lower.tail = FALSE
                )
  freq_granges <- lowcountgranges[res < binom.p]
  freqrejects <- length(lowcountgranges) - length(freq_granges)
  if(append == TRUE) {
    output_granges <- freq_granges
  }
  ##loginfo("...done")
############################################################################
  if (append != TRUE) {
    ##Performing the Fisher's Test Filters
    if(fet_2x2 == TRUE){
      message("Filtering granges by Fishtest pvals...")
      two_by_2_pvals <- list()
      temp_mat1 <- cbind(values(freq_granges)$count.pos.ref,
                         values(freq_granges)$count.pos,
                         values(freq_granges)$count.neg.ref,
                         values(freq_granges)$count.neg)
      x1 <- vector(mode = 'numeric', length = dim(temp_mat1)[1])
      for(j in 1:dim(temp_mat1)[1]) {
        if (sum(as.numeric(temp_mat1[j,])) < 90000) {
          x1[j] <- fisher.test(matrix(temp_mat1[j,], ncol = 2, byrow = F), workspace = 1000000)$p.val
        } else {
          x1[j] <- fisher.test(matrix(temp_mat1[j,], ncol = 2, byrow = F), workspace = 1000000, simulate.p.value = TRUE)$p.val
        }
      }
      index <- which(x1 >= pval)
      two_by_two_filter_gr <- freq_granges[index]
      fishers_2by2_rejects <- length(freq_granges) - length(two_by_two_filter_gr)
    }else{
        two_by_two_filter_gr <- freq_granges
    }
    if(fet_3x2 == TRUE){
      warning("Hard-coded to expect tally2GR was given three bins. Generalize this")
      three_by_2_pvals <- list()
      temp_mat <- cbind(
                        values(two_by_two_filter_gr)[,16],
                        values(two_by_two_filter_gr)[,13],
                        values(two_by_two_filter_gr)[,17],
                        values(two_by_two_filter_gr)[,14],
                        values(two_by_two_filter_gr)[,18],
                        values(two_by_two_filter_gr)[,15]
                        )
      x <- vector(mode = 'numeric', length = dim(temp_mat)[1])
      for(i in 1:dim(temp_mat)[1]) {
        if (sum(as.numeric(temp_mat[i,])) < 90000) {
          x[i] <- fisher.test(matrix(temp_mat[i,], ncol = 3, byrow = F), workspace = 1000000)$p.val
        } else {
          x[i] <- fisher.test(matrix(temp_mat[i,], ncol = 3, byrow = F), workspace = 1000000, simulate.p.value = TRUE)$p.val
        }
      }
      index2 <- which(x >= pval)
      output_granges <- two_by_two_filter_gr[index2]
      fishers_3by2_rejects <- length(two_by_two_filter_gr) - length(output_granges)
    }else{
      output_granges <- two_by_two_filter_gr
    }
  }
  ##loginfo("...done")
  message("outputing result...")
  output <- list()
  if (append) {
    output[[1]] <- granges1
  } else {
    output[[1]] <- granges
  }
  output[[2]] <- output_granges
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
  names(output) <- c("raw_granges", "filtered_granges", "N_nucleotide_rejects", "low_count_rejects", "frequency_rejects", "fisher_test_rejects_2by2", "fisher_test_rejects_3by2")
  return(output)
  #saveWithID(output, "raw_and_filtered_granges", prepend_str, save_dir=file.path(save_dir, "Breaks_RData"), compress=FALSE)
  ##loginfo("...done")
  message("done!")
}

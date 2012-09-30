
 transparent <- function (cols, alpha = 0.4) {
    col_mat <- col2rgb(cols, alpha = TRUE)
    col_mat["alpha", ] <- alpha * 255
    result <- rgb(col_mat["red", ],
                  col_mat["green", ],
                  col_mat["blue",],
                  col_mat["alpha", ], maxColorValue = 255)
    names(result) <- names(cols)
    return(result)
  }

p.L <- function(x, n, alpha) {
  qbeta(alpha, x, n - x + 1)
}
p.U <- function(x, n, alpha) {
  qbeta(1 - alpha, x + 1, n - x)
}


.plotVarChrom <- function(normal_gr, tumor_gr, normal_cov_gr, chrom, freq = TRUE){
  overN <- findOverlaps(tumor_gr, normal_cov_gr)
  cat(length(tumor_gr[-queryHits(overN)]), "\n", sep = "")
  tumor_gr <- tumor_gr[queryHits(overN)]
  tumor_gr <- tumor_gr[(values(tumor_gr)$count)>3]
  over1 <- which(paste(values(tumor_gr)$location,values(tumor_gr)$alt)
                 %in%
                 paste(values(normal_gr)$location, values(normal_gr)$alt))
  tumor_and_norm <- tumor_gr[over1]
  tumor_specific <- tumor_gr[-over1]
  chrom_gr <- tumor_and_norm[seqnames(tumor_and_norm)==chrom]
  
  par(mfrow = c(2,1))
  plot(start(chrom_gr)/10e5,
       values(chrom_gr)$count/(values(chrom_gr)$count.total), pch = 19,
       cex = .5, col = "red", ylim = c(0,1),
       xlim = c(0,max(start(chrom_gr)/10e5)),
       xlab = paste(chrom, " position in Mb"),
       ylab = "Read frequency of Alt allele")
  
  chrom_ts_gr <- tumor_specific[seqnames(tumor_specific)==chrom]
  
  par(new=T)
  plot(start(chrom_ts_gr)/10e5,
       values(chrom_ts_gr)$count/(values(chrom_ts_gr)$count.total),
       pch = 19, cex = 1, ylim = c(0,1),
       xlim = c(0,max(start(chrom_gr)/10e5)), axes =F, xlab=NA, ylab =NA)

  if(freq==TRUE){
    c1 <- hist(values(chrom_ts_gr)$count/(values(chrom_ts_gr)$count.total), breaks = seq(0, 1, by = 0.05), plot=FALSE)$count
    c2 <- hist(values(chrom_gr)$count/(values(chrom_gr)$count.total), breaks = seq(0, 1, by = 0.05), plot=FALSE)$count
     barplot(rbind(c1, c2), beside=T,
            main = "Tumor-specific variant frequency",
            ylab = "Fraction", xlab = "Total count",
             legend = c("Tumor specific", "Tumor and normal"))
  }else{
    
    c1 <- hist(values(chrom_gr)$count.total, breaks = seq(0, max(values(chrom_gr)$count.total)+100, by=100), plot = FALSE)$count
    c2<- hist(values(chrom_ts_gr)$count.total, breaks = seq(0,  max(values(chrom_gr)$count.total)+100, by=100), plot = FALSE)$count
    barplot(rbind(c1/sum(c1), c2/sum(c2)), beside=T,
            main = "Tumor-specific variant coverage",
            ylab = "Count", xlab = "Total count")
  }
}

 plotVarGenome <- function(normal_gr, normal_raw, tumor_gr, normal_cov,
                           tumor_rna=NULL, freq =FALSE, SRC = NULL,
                           main = "NULL", expressed = FALSE, mosaic = TRUE,
                           bycount = FALSE){
   
   tumor_gr <- tumor_gr[(values(tumor_gr)$count)>=2]
   message("finding normal failed vs normal low freq")
   over1 <- which(paste(values(normal_raw)$location,values(normal_raw)$alt)
                  %in%
                  paste(values(normal_gr)$location, values(normal_gr)$alt))
   
   not_called <- normal_raw[-over1]
   
   normal_low <-not_called[values(not_called)$count/values(not_called)$count.total <.04023]
   
   failed_normal <- not_called[values(not_called)$count/values(not_called)$count.total >=.04023]
   
   ## make normal to test
   
   normal <- c(normal_gr, failed_normal)
   
   message( "find the variants called in tumor and in normal")
   over1 <- which(paste(values(tumor_gr)$location,values(tumor_gr)$alt)
                  %in%
                  paste(values(normal_gr)$location, values(normal_gr)$alt))
   
   tumor_and_normal <- tumor_gr[over1]
   
   message("find the tumor specific calls")
   over1 <- which(paste(values(tumor_gr)$location,values(tumor_gr)$alt)
                  %in%
                  paste(values(normal_raw)$location, values(normal_raw)$alt))
   
   tumor_specific <- tumor_gr[-over1]
   
   over2 <- which(paste(values(tumor_gr)$location,values(tumor_gr)$alt)
               %in%
                  paste(values(normal_low)$location, values(normal_low)$alt))
   
   over3 <- which(paste(values(normal_low)$location,values(normal_low)$alt)
                  %in%
                  paste(values(tumor_gr)$location, values(tumor_gr)$alt))
   
   low_in_tumor <- normal_low[over3]
   
   tumor_in_low <- tumor_gr[over2]
   
   grl <- split(tumor_specific, as.factor(seqnames(tumor_specific)))
#### Need to extend the coverage vectors to the max position of the grl for each chromosome
   
   ends <- lapply(grl, function(x) { max(end(x))})
   
   cov <- lapply(names(ends), function(x) {
    if(ends[[x]] > length(normal_cov[[x]])){
      extend = Rle(values=0, length= ends[[x]]- length(normal_cov[[x]]))
      normal_cov[[x]] = c(normal_cov[[x]], extend)
    }
    return(normal_cov[[x]])
  })         
  names(cov) <- names(grl)
  tumor_specific<- unlist(grl)
  p.lower <-  qbeta(0.125,
                    values(tumor_specific)$count,
                    values(tumor_specific)$count.total -
                    values(tumor_specific)$count +1)
  ##p.lower <- values(tumor_specific)$count/values(tumor_specific)$count.total
  names <-as.list(names(grl))
  n <- lapply(names, function(x) {as.numeric(cov[[x]][ranges(grl[[x]])])})
  n <- unlist(n)
  #p.lower2 <- values(tumor_specific)$count/values(tumor_specific)$count.total
   p.bin <- pbinom(0, n, p.lower)
   cat("Insuficient  power in normal ", sum(p.bin>.01), "\n", sep = "")
   tumor_specific<-tumor_specific[p.bin<0.01]
   p.lower <-  qbeta(0.125,
                     values(tumor_in_low)$count,
                     values(tumor_in_low)$count.total -
                     values(tumor_in_low)$count +1)
   p.bin <- pbinom(values(low_in_tumor)$count, values(low_in_tumor)$count.total, p.lower)
   add <- tumor_in_low[p.bin<0.01]
   
   tumor_specific <- c(tumor_specific, add)

  if(expressed==TRUE){
    over1 <- which(paste(values(tumor_specific)$location,values(tumor_specific)$alt)
               %in%
               paste(values(tumor_rna)$location, values(tumor_rna)$alt))
    cat("Unexpressed tumor variants ", length(tumor_specific)-length(over1), "\n", sep ="")
    tumor_specific <- tumor_specific[over1]
  }
   p<- values(tumor_and_normal)$count/(values(tumor_and_normal)$count.total)
   ##x <- values(tumor_and_normal)$count
   ##cov <- values(tumor_and_normal)$count.total
   ##ci = p.U(x, cov, 0.025) -p.L(x, cov, 0.025) 
   if(mosaic ==TRUE){
     nf <- layout(matrix(c(1,1,2,3), 2,2, byrow=TRUE), widths = c(1.5, 1.5), respect=TRUE)
   }
   if(freq ==TRUE){
      nf <- layout(matrix(c(1,1,2,2), 2, 2, byrow=TRUE), respect=TRUE)
    }
   chroms <- cumsum(as.numeric(seqlengths(Hsapiens)[1:24]))
   chroms[2:24] <- chroms[1:23]
   chroms[1] <- 0
   names(chroms) <- names(seqlengths(Hsapiens)[1:24])
   values(tumor_and_normal)$genome_pos <- start(tumor_and_normal)+chroms[as.character(seqnames(tumor_and_normal))]
   values(tumor_specific)$genome_pos <- start(tumor_specific)+chroms[as.character(seqnames(tumor_specific))]
  
   plot(values(tumor_and_normal)$genome_pos/10e5,
        p, pch = 19,
        cex = .1, col = "dark grey", ylim = c(0,1),
        xlim = c(0,max(values(tumor_and_normal)$genome_pos/10e5, na.rm=T)),
        xlab = paste("Position on genome"),
        ylab = "Read frequency of Alt allele",
        main = SRC)
   
   
   par(new=T)
   plot(values(tumor_specific)$genome_pos/10e5,
        values(tumor_specific)$count/(values(tumor_specific)$count.total),
        pch = 19, cex = .4, ylim = c(0,1),
        xlim = c(0,max(values(tumor_and_normal)$genome_pos/10e5, na.rm =T)), axes =F, xlab=NA, ylab =NA, col="orangered")
   abline(v = chroms/10e5)
   if(mosaic == TRUE){
     mosaicplot(table(as.character(values(tumor_and_normal)$ref), as.character(values(tumor_and_normal)$alt)),
                col = c("firebrick3", "dodgerblue3", "forestgreen", "goldenrod3"),
                cex.axis = 2,
                main = "Germline variants")
     mosaicplot(table(as.character(values(tumor_specific)$ref), as.character(values(tumor_specific)$alt)),
                col = c("firebrick3", "dodgerblue3", "forestgreen", "goldenrod3"),
                cex.axis = 2,
                main = "Tumor specific variants")
   }
   if(freq==TRUE){
     c1 <- hist(values(tumor_specific)$count/(values(tumor_specific)$count.total), breaks = seq(0, 1, by = 0.01), plot=FALSE)$count
     c2 <- hist(values(tumor_and_normal)$count/(values(tumor_and_normal)$count.total), breaks = seq(0, 1, by = 0.01), plot=FALSE)$count
     barplot(rbind(c1/sum(c1), c2/sum(c2)), beside=T,
             main = "Variant frequencies",
             ylab = "Fraction", xlab = "Alternate read frequency",
             legend = c("Tumor specific", "Tumor and normal"),
             names =seq(.01, 1, by = 0.01),
             las=2)
   }
   if(bycount ==TRUE){
     c1 <- hist(values(tumor_and_normal)$count.total, breaks = seq(0, max(values(tumor_and_normal)$count.total)+100, by=100), plot = F)$count
     c2<- hist(values(tumor_specific)$count.total, breaks = seq(0,  max(values(tumor_and_normal)$count.total)+100, by=100), plot = F)$count
     barplot(rbind(c1/sum(c1), c2/sum(c2)), beside=T,
             main = main,
             ylab = "Count", xlab = "Total count")
   }
 }


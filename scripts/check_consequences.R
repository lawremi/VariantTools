x <- Lgr[["SAM587376"]]

sams <- names(kras)
sams_list <- as.list(sams)
with_cons <- lapply(sams_list, function(x){
  gr <- kras[[x]]
  if(file.exists(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/TN/"
                       , x,"_variant_annotations.txt", sep =""))){
    tab <- read.table(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/TN/"
                            , x,"_variant_annotations.txt", sep =""), sep = "\t")
  }else{
    return(gr)
  }
  if(length(gr)>0){
    values(gr)$cons <- NA
    pos_x <- paste(sub("chr", "", seqnames(gr)), start(gr), sep = ":")
    for(i in 1:length(pos_x)){
      cons <- unique(as.character(tab[which(tab$V2 == pos_x[i]), 7]))
      if(length(cons) ==1){
        values(gr[i])$cons<- cons
      }else{
        cons <- paste(cons, collapse=",")
        values(gr[i])$cons<- cons
      }
    }
    return(gr)
  }else{
    return(gr)
  }
})

final <- lapply(sams_list, function(sam){
  x <- with_cons[[sam]]
  if(length(x)>0){
    t <- cbind(sam, as.character(seqnames(x)),
               start(x), as.character(values(x)$ref),
               as.character(values(x)$read),values(x)$count,
               values(x)$count.total,values(x)$cons)
    return(t)
  }else{
    return(NULL)
  }
})


merge_var <- function(gr_list){
  res <- list()
  list_breaks <- PartitioningByWidth(elementLengths(gr_list))
  ul_result <- unlist(gr_list)
  names(ul_result)<-NULL
  res$raw_granges<-do.call(c,ul_result[start(list_breaks)])
  res$filtered_granges<-do.call(c,ul_result[start(list_breaks)+1])
  res$N_nucleotide_rejects <- do.call(sum,ul_result[start(list_breaks)+2])
  res$low_count_rejects <- do.call(sum,ul_result[start(list_breaks)+3])
  res$frequency_rejects <- do.call(sum,ul_result[start(list_breaks)+4])
  res$fisher_test_rejects_2by2 <- do.call(sum,ul_result[start(list_breaks)+5])
  res$read_pos_rejects <- do.call(sum,ul_result[start(list_breaks)+6])
  return(res)
}

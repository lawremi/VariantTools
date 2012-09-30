function(Qscore, count, offset){
  pb <- 1-(10^(-(Qscore-offset)/10))
  prob_count <- sum(pb*count)
  return(prob_count)
}

callWildtype <- function(cov, var_gr, power=0.999, plower=0.2, perror=1/1000) {
  callable <- callCallable(cov, plower = plower, perror = perror, power = power)
  wildtype <- callable
  wildtype[!callable] <- NA
  rl <- split(ranges(var_gr), factor(seqnames(var_gr), names(cov)))
  wildtype[rl] <- FALSE
  wildtype
}

byCov <- function(sam){
  som <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP2.0_var_rerun/TN/",
                  sam , "/",
                  sam,
                  ".sample_specific_variants_granges.RData", sep = "")))
  som <- som[values(som)$count/values(som)$count.total >.1]
  cov <- get(load(paste("/gne/research/data/cgp/2.0/", sam, "/exome_merged/RData/",sam,".cov.RData", sep = "")))
  X <- slice(cov, lower = 8, upper =12)
  XX <- slice(cov, lower = 18, upper = 22)
  XXX <- slice(cov, lower = 28, upper = 32)
  gr10x <- as(X, "GRanges")
  gr20x <- as(XX, "GRanges")
  gr30x <- as(XXX, "GRanges")
  over <- findOverlaps(som, gr10x)
  in10 <- som[queryHits(over)]
  over <- findOverlaps(som, gr20x)
  in20 <- som[queryHits(over)]
  over <- findOverlaps(som, gr30x)
  in30 <- som[queryHits(over)]
  x <- c(length(in10)/(sum(width(gr10x))/10e5), length(in20)/(sum(width(gr20x))/10e5), length(in30)/(sum(width(gr30x))/10e5))
  barplot(x, main = sam, ylab = "Tumor-specific mutations per Mb", names = c("10X", "20X", "30X"))
}

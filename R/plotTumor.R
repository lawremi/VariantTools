## plotVariantsAlongGenome <- function(x, ...) {
##   plotGrandLinear(x, aes(y = count / count.total), ...,
##                   ylab = "Read frequency of ALT allele")
## }

## plotTNVariantsAlongGenome <- function(tumor, normal, ...) {
##   variants <- mstack(tumor = tumor, normal = normal, .indName = "tumor")
##   plotVariantsAlongGenome(variants, aes(color = tumor, size = tumor), ...) +
##     scale_color_manual(values = c("darkgray", "orangered")) +
##       scale_size_discrete(range = c(0.1, 0.4))
## }

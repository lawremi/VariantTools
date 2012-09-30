plotTiTv <- function(gr, main = "TiTv"){
  titv <- table(as.character(values(gr)$ref), as.character(values(gr)$alt))
  col = c("firebrick3", "dodgerblue3", "forestgreen", "goldenrod3")
  mosaicplot(titv, col = col, main = main)
}

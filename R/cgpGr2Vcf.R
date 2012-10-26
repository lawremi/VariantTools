## Convert Genentech Variant GRanges into a VCF 

## How we will populate a VCF object:

## rowData: just the ranges
## colData: N/A
## exptData: 'header' element of class VCFHeader, with 'header' slot
##           being a DataFrameList, with a FORMAT element, where rownames => ID,
##           and the other columns are Number, Type and Description.
## fixed: REF and ALT
## info: N/A
## geno: list of single-column matrices for AR, RR, DP, AAP, RAP

variantGR2Vcf <- function(x, sample.id, project = NULL) {
  location <- factor(x$location, unique(x$location))
  locationRle <- Rle(location)
  
  rowData <- x[start(locationRle)]
  mcols(rowData) <- NULL

  colData <- DataFrame(Samples = 1L, row.names = sample.id)
  
  meta_vector <- c(source = paste("VariantTools",
                     packageVersion("VariantTools")),
                   phasing = "unphased", fileformat = "VCFv4.1",
                   project = project)
  meta_header <- DataFrame(Value = meta_vector, row.names = names(meta_vector))
  format_header <- read.csv(system.file("vcf", "header.csv",
                                        package = "VariantTools"))
  header <- VCFHeader(reference = seqlevels(x), samples = sample.id,
                      header = DataFrameList(META = meta_header,
                        FORMAT = as(format_header, "DataFrame")))
  exptData <- SimpleList(header = header)

  fixed <- DataFrame(REF = DNAStringSet(x$ref[start(locationRle)]),
                     ALT = split(DNAStringSet(x$alt), location),
                     QUAL = rep.int(NA_real_, length(rowData)),
                     FILTER = rep.int(NA_character_, length(rowData)))

  genoMatrix <- function(v) {
    matrix(v, nrow = nrow(fixed), 1)
  }
  
  geno <-
    SimpleList(AR= genoMatrix(split(x$count, location)),
               RR= genoMatrix(x$count.ref[start(locationRle)]),
               DP= genoMatrix(x$count.total[start(locationRle)]),
               AAP=genoMatrix(split(as.integer(x$count > 0L), location)),
               RAP=genoMatrix(as.integer(x$count.ref[start(locationRle)] > 0L)))

  VCF(rowData = rowData, colData = colData, exptData = exptData, fixed = fixed,
      geno = geno)
}

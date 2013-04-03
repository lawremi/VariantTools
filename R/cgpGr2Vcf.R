## Convert Genentech Variant GRanges into a VCF 

variantGR2Vcf <- function(x, sample.id, project = NULL,
                          genome = unique(GenomicRanges::genome(x)))
{
  if (missing(sample.id) || !isSingleString(sample.id))
    stop("'sample.id' must be provided as a single, non-NA string")
  if (!is.null(project) && !isSingleString(project))
    stop("'project', if provided, must be a single, non-NA string")
  
  rowData <- x
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

  x <- normalizeIndelAlleles(x, genome)
  
  fixed <- DataFrame(REF = DNAStringSet(x$ref),
                     ALT = as(as.character(x$alt), "List"),
                     QUAL = rep.int(NA_real_, length(rowData)),
                     FILTER = rep.int(NA_character_, length(rowData)))

  genoMatrix <- function(v) {
    matrix(v, nrow = nrow(fixed), 1)
  }

  posFactor <- rep(seq_len(length(x)), 2)
  alleleDepth <- c(x$count.ref, x$count)
  allelePresent <- as.integer(c(x$count.ref > 0, x$count > 0L))

  geno <-
    SimpleList(AD = genoMatrix(split(alleleDepth, posFactor)),
               DP = genoMatrix(x$count.total),
               AP = genoMatrix(split(as.integer(allelePresent), posFactor)))

  VCF(rowData = rowData, colData = colData, exptData = exptData, fixed = fixed,
      geno = geno)
}

normArgGenome <- function(x) {
  if (isSingleString(x))
    x <- GmapGenome(x)
  else if (!is(x, "GmapGenome"))
    stop("'genome' must be either a 'GmapGenome' object or ",
         "a single string identifying one")
  x
}

normalizeIndelAlleles <- function(x, genome = unique(GenomicRanges::genome(x))) {
  is.indel <- nchar(x$ref) == 0L | nchar(x$alt) == 0L
  if (any(is.indel)) {
    genome <- normArgGenome(genome)
    indels <- x[is.indel]
    indels <- shift(indels, -1)
    anchor <- getSeq(genome, indels)
    indels$ref <- paste0(anchor, indels$ref)
    indels$alt <- paste0(anchor, indels$alt)
    x[is.indel] <- indels
  }
  x
}

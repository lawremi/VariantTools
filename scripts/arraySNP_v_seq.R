library(genoset)
library(BigMatrix)
library(CGPtools)
library(VariantTools)
library(BSgenome.Hsapiens.UCSC.hg19)
source("~/svn/Packages/VariantTools/scripts/getAlleleByLocus_mod.R")

meta<-getCGPMetadata(JSON="http://research.gene.com/cgp/released.json")
## ds = readGenoSet("/gne/research/data/cgp/SNPArray.V4/illumina_2.5M.genoset.genotypes.and.raw.RData")
ds = readCGPSNPArray(version="Genotype")   # Very recent addition
my.levels = levels(ds[ , , "geno"])

het.levels = c("AC","AG","AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG") # or something

locs_gr <- as(locData(ds), "GRanges")
seqlevels(locs_gr) <- paste("chr", seqlevels(locs_gr), sep = "")  
locs_gr <- locs_gr[!duplicated(paste(seqnames(locs_gr), start(locs_gr), sep = ""))]

plot_SNPS("SAM587284")

plot_SNPS<- function(SAM){
locs_gr <- as(locData(ds), "GRanges")
seqlevels(locs_gr) <- paste("chr", seqlevels(locs_gr), sep = "")
cov <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS121/", SAM, "/results/",SAM,".coverage.RData", sep = "")))

gr1 <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/",
                      SAM,
                      "/results/",
                      SAM,
                      ".variants_granges.RData", sep = "")))
seqlevels(gr1) <- paste("chr", seqlevels(gr1), sep = "")
names(cov) <- paste("chr", names(cov), sep = "")
names <- seqlevels(locs_gr)
grl <- split(locs_gr, as.factor(seqnames(locs_gr)))
n <- lapply(names, function(x) { as.numeric(cov[[x]][ranges(grl[[x]])]) })
values(locs_gr)$cov <- unlist(n)
int <- which(pData(ds)$SAMPLE_ID ==sub("SAM", "", SAM))
geno <- ds[ , int, "geno"]
values(locs_gr)$geno <-geno
values(locs_gr)$count <-0
over <- findOverlaps(gr1, locs_gr)
values(locs_gr[subjectHits(over)])$count <- values(gr1[queryHits(over)])$count
values(locs_gr)$alt <- NA
values(locs_gr[subjectHits(over)])$alt <- as.character(values(gr1[queryHits(over)])$read)
seq <- getSeq(Hsapiens, locs_gr)
values(locs_gr)$ref<-as.character(seq)
##values(locs_gr[subjectHits(over)])$ref <- as.character(values(gr1[queryHits(over)])$ref)

values(locs_gr)$isHet <- as.character(values(locs_gr)$geno)%in% het.levels
values(locs_gr)$homRef <- (values(locs_gr)$geno ==paste(values(locs_gr)$ref, values(locs_gr)$ref, sep=""))
gr2 <- get(load(paste("/gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/",
                      SAM,
                       "/results/",
                       SAM,
                       ".filtered_variants_granges.RData", sep = "")))
seqlevels(gr2) <- paste("chr", seqlevels(gr2), sep = "")

with_cov <- locs_gr[values(locs_gr)$cov>15]
SNPS <- file.path(system.file(package = "VariantTools", mustWork=TRUE),"extdata/exon_snps_ill25.RData")
SNPs <- get(load(SNPS))
over <- findOverlaps(SNPs, with_cov)
with_cov[subjectHits(over)]->ex

het <- ex[values(ex)$isHet]
hom <- ex[!values(ex)$isHet]
hom_ref <- hom[which(values(hom)$homRef)]
hom_alt <-hom[which(!values(hom)$homRef)]

over1 <- findOverlaps(het, gr2)
called_het <- het[queryHits(over1)]
not_called_het <- het[-queryHits(over1)]
over2 <- findOverlaps(hom_ref, gr2)
called_hom_ref <- hom_ref[queryHits(over2)]
not_called_hom_ref <- hom_ref[-queryHits(over2)]
over3 <- findOverlaps(hom_alt, gr2)
called_hom_alt <- hom_alt[queryHits(over3)]
not_called_hom_alt <- hom_alt[-queryHits(over3)] 

par(mfrow = c(1,3))
plot(log10(values(het)$count), log10(values(het)$cov), xlab = "log10 alt read count", ylab = "log10 total count", main = "Heterozygous SNPs by array", pch = NA)
points(log10(values(called_het)$count), log10(values(called_het)$cov), col = "blue", pch = 19)
points(log10(values(not_called_het)$count), log10(values(not_called_het)$cov), col = "red", pch = 19, cex = 0.4)
abline(0,1)
abline(.301,1)
plot(log10(values(hom_alt)$count), log10(values(hom_alt)$cov), xlab = "log10 alt read count", ylab = "log10 total count", main = "Homozygous Alt SNPs by array", pch = NA)
points(log10(values(called_hom_alt)$count), log10(values(called_hom_alt)$cov), col = "blue", pch = 19)
points(log10(values(not_called_hom_alt)$count), log10(values(not_called_hom_alt)$cov), col = "red", pch = 19, cex = 0.4)
abline(0,1)
abline(.301,1) 
plot(log10(values(hom_ref)$count), log10(values(hom_ref)$cov), xlab = "log10 alt read count", ylab = "log10 total count", main = "Homozygous Ref SNPs by array", pch = NA)
points(log10(values(called_hom_ref)$count), log10(values(called_hom_ref)$cov), col = "blue", pch = 19)
points(log10(values(not_called_hom_ref)$count), log10(values(not_called_hom_ref)$cov), col = "red", pch = 19, cex = 0.4)
abline(0,1)
abline(.301,1)
legend(x=2.5, y = 3.5, col = c("blue", "red"), legend=c("Called SNV", "No call"), pch = 19)
}

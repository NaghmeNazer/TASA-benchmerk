library(data.table)
library(GEOquery)
library(minfi)
library(valr)
library(doParallel)

## function to find DMRs with minfi bumphunter method.
minfi_bumphunter <- function(B.mat, Pheno, cutoff = 0.3, maxGap = 50, permutation_number = 1000, cores = 15){
  GRset <- makeGenomicRatioSetFromMatrix(B.mat, what = "Beta")
  designMatrix <- model.matrix(~ Pheno)
  totalCores = detectCores()
  cluster <- makeCluster(cores) 
  registerDoParallel(cluster)
  dmrs <- bumphunter(GRset, design = designMatrix, cutoff = cutoff, maxGap = maxGap, B=permutation_number, type="Beta")
  return(dmrs[["table"]])
}

## read beta values matrix
data <- fread("../GSE87053/GSE87053_b.csv", data.table = F, header = T)
rownames(data) <- data$V1
data$V1 <- NULL
data.mat <- as.matrix(data)
rownames(data.mat) <- rownames(data)

## gets the phenotype of samples
gse <- getGEO(filename="../GSE87053/GSE87053_series_matrix.txt.gz", getGPL = F)
pheno_data <- gse@phenoData@data[,c("geo_accession","source_name_ch1")]
pheno <- pheno_data[colnames(data), "source_name_ch1"]
table(pheno)

# run bumphunter
bump_dmrs <- minfi_bumphunter(data.mat, pheno)
fwrite(bump_dmrs, "../GSE87053/bump_dmrs.csv")

#retain only significant DMRs with p-val < 0.05
bump_dmrs.sig <- bump_dmrs[bump_dmrs$p.value < 0.05,]

### read paper candidate genes
genes <- fread("../GSE87053/confirmed_genes.csv", data.table = F)
genes_bed <- genes[,c("chrom", "start", "end")]

# make bumphunter bed file
bump_dmrs_bed <- bump_dmrs.sig[,c("chr", "start","end")]
names(bump_dmrs_bed)[1] <- "chrom"
bump_dmrs_bed$chrom <- as.numeric(gsub("chr", "", bump_dmrs_bed$chrom))
bump_dmrs_bed <- bump_dmrs_bed[!is.na(bump_dmrs_bed$chrom),]

# intersect bumphunter results with confirmed genes
bump_genes_interval <- bed_intersect(genes_bed, bump_dmrs_bed)
length(unique(bump_genes_interval$start.x))
bump_genes_not_interval <- bed_intersect(genes_bed, bump_dmrs_bed, invert = T)
length(unique(bump_genes_not_interval$start))

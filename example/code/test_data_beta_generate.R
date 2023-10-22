### Load Packages
library(GEOquery)
library(minfi)
library(data.table)
library(ggplot2)

### Functions
## function to remove samples with low medial intensity in M and U signals
sample.qc <- function(M, U, thresh = 10){
  Mmed <- log2(apply(M,2,median))
  Umed <- log2(apply(U,2,median))
  !(Mmed < 10 & Umed <10)
  good.samples <- !(Mmed < 10 & Umed <10)
  return(good.samples)
}

## function to generate beta values from M and U signals
intensities.to.beta <- function(M, U, alpha = 100){
  B <- M/(M+U+alpha)
  return(B)
}

## function to remove probes with high detection p value
bad.detected.probes <- function(P){
  if (is.null(P)) return(NULL)
  else return(unique(rownames(which(P > 0.05, arr.ind = T))))
}

## function to remove probes with negative intensity value
negative.probes <- function(M, U){
  return(union(unique(rownames(which(M < 0, arr.ind = T))), unique(rownames(which(U < 0, arr.ind = T)))))
}

## function to remove probes containing a know SNP with high allele freq
rs_probes <- function(){
  rs <- read.csv('../info_datasets/rs_af0.05.csv')
  return(unique(rs$X))
}

## function to remove probes mapping to multiple loci in genome
non_specific.probes <- function(){
  ps <- read.csv('../info_datasets/non_specific_sites.csv')
  return(unique(ps[,2]))
}

## function to remove probes with constant value in all samples
constant_value.probes <- function(M, U){
  ps <- c()
  row.sds <- rowSds(M)
  ps <- c(ps, rownames(M)[which(row.sds == 0)])
  row.sds <- rowSds(U)
  ps <- c(ps, rownames(M)[which(row.sds == 0)])
  ps <- unique(ps)
  return(ps)
}

## function to aggregate all bad quality probes
deleted.probes <- function(M,U,P = NULL){
  high.p.probes <- bad.detected.probes(P)
  neg.probes <- negative.probes(M, U)
  constant.val.probes <- constant_value.probes(M, U)
  nspec.probes <- non_specific.probes()
  snp.probes <- rs_probes()
  deleted.probes <- unique(c(high.p.probes, neg.probes, 
                             constant.val.probes, nspec.probes, snp.probes))
  return(deleted.probes)
}

## function to preprocess and generate the beta value based on M, U and P matrices
preprocess_to_beta <- function(M, U, P = NULL){
  bad.probes <- deleted.probes(M,U,P)
  good.samples <- sample.qc(M, U)
  B <- intensities.to.beta(M,U)
  B <- B[!(rownames(B) %in% bad.probes),]
  B <- B[,good.samples]
  return(B)
}

### Main
### Read GSE87053
#decompress idats
untar("../GSE87053/GSE87053_RAW.tar", exdir = "../GSE87053/idat")
#list files
head(list.files("../GSE87053/idat", pattern = "idat"))
idatFiles <- list.files("../GSE87053/idat", pattern = "idat.gz$", full = TRUE)
#decompress individual idat files
sapply(idatFiles, gunzip, overwrite = TRUE)
#read idats and create RGSet and MSet
RGSet <- read.metharray.exp("../GSE87053/idat")
MSet <- preprocessRaw(RGSet)
# get M and U signal intensities
M <- getMeth(MSet)
U <- getUnmeth(MSet)
# preprocess and generate beta value matrix
B <- preprocess_to_beta(M, U)
colnames(B) <- gsub("_.*", "",colnames(B))
# write the beta value matrix for furthet analysis
fwrite(data.frame(B), "../GSE87053/GSE87053_b.csv", row.names = T, col.names = T)

# Perform PCA
# read series matrix file to get the pheno informtaion
gse <- getGEO(filename="../GSE87053/GSE87053_series_matrix.txt.gz", getGPL = F)
pheno_data <- gse@phenoData@data[,c("geo_accession","source_name_ch1")]
pheno <- pheno_data[colnames(B), "source_name_ch1"]
table(pheno)
# perform PCA and plot the first two PCs
plot_df <- data.frame(t(B))
pca_res <- prcomp(plot_df, scale = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])
head(pca_res.x)
pdf("../GSE87053/GSE87053-qced-PCA-plot.pdf")
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = pheno)) + geom_point() + ggtitle("GSE87053 Dataset") + theme_bw()
dev.off()

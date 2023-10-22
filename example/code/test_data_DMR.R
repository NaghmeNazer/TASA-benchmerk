library(data.table)
library(GEOquery)
library(ENmix)
library(valr)
library(doParallel)

## read illumina 450K manifest file
manifest <- fread("../info_datasets/450k_manifest.csv", data.table = F)
rownames(manifest) <- manifest$IlmnID 

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

##read DMRcate DMPs
dmps <- fread("../GSE87053/dmrcate_dmps.csv")
names(dmps)[1] <- "probe"
names(dmps)[4] <- "p"
dmps$chr <- manifest[dmps$probe, "CHR"]
dmps$start <- manifest[dmps$probe, "MAPINFO"]
dmps$end <- dmps$start + 1
dmps <- dmps[,c("chr", "start", "end", "p", "probe")]

# run ipdmr
ipdmr(data=dmps,seed=0.05, region_plot = F, mht_plot = F)
ipdmr_dmrs <- fread("resu_ipdmr.csv", data.table = F)
fwrite(ipdmr_dmrs, "../GSE87053/ipdmr_dmrs.csv")
gc()
ipdmr_dmrs_morethanone <- ipdmr_dmrs[ipdmr_dmrs$nprobe > 1,]

### read paper candidate genes
genes <- fread("../GSE87053/confirmed_genes.csv", data.table = F)
genes_bed <- genes[,c("chrom", "start", "end")]

# make ipdmr bed file
ipdmr_dmrs_bed <- ipdmr_dmrs[,c("chr", "start", "end")]
names(ipdmr_dmrs_bed)[1] <- "chrom"
ipdmr_dmrs_bed$chrom <- as.numeric(ipdmr_dmrs_bed$chrom)
ipdmr_dmrs_bed <- ipdmr_dmrs_bed[!is.na(ipdmr_dmrs_bed$chrom),]

# intersect ipdmr results with confirmed genes
ipdmr_genes_interval <- bed_intersect(genes_bed, ipdmr_dmrs_bed)
length(unique(ipdmr_genes_interval$start.x))
ipdmr_genes_not_interval <- bed_intersect(genes_bed, ipdmr_dmrs_bed, invert = T)
length(unique(ipdmr_genes_not_interval$start))

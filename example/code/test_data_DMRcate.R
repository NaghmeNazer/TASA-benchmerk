library(DMRcate)
library(data.table)
library(GEOquery)

## read beta value matrix
data <- fread("../GSE87053/GSE87053_b.csv", data.table = F, header = T)
rownames(data) <- data$V1
data$V1 <- NULL
data.mat <- as.matrix(data)
rownames(data.mat) <- rownames(data)

## identify the phenotype of each sample
gse <- getGEO(filename="../GSE87053/GSE87053_series_matrix.txt.gz", getGPL = F)
pheno_data <- gse@phenoData@data[,c("geo_accession","source_name_ch1")]
pheno <- pheno_data[colnames(data), "source_name_ch1"]
table(pheno)
design <- model.matrix(~factor(pheno))

## run DMRcate to get the dmps and save it.
dmps <- cpg.annotate(datatype = "array", data.mat, what = "Beta", arraytype = "450K", 
                     analysis.type = "differential", design, coef = 2)
dmps.df <- cbind(dmps@ranges@ranges@NAMES,data.frame(dmps@ranges@elementMetadata@listData))
dmps.df.sig <- dmps.df[dmps.df$is.sig,]
fwrite(dmps.df.sig, "../GSE87053/dmrcate_dmps.csv")

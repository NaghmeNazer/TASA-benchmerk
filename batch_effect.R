## Simulation datasets from different sources (GSE120610, GSE131989, GSE134429, GSE184269, GSE56046)
## are loaded and batch effect removed with combat from sva package.
## batch effect cirrected/not corrected datasets are saved for downstream analysis
library(data.table)
library(ggplot2)
library(sva)

## load different datasets
monos_1 <- fread('/home/sepehri/simulation/Mono/beta_chr1_GSE120610.csv', data.table = F)
rownames(monos_1) <- monos_1$V1
monos_1$V1 <- NULL

cd8s_1 <- fread('/home/sepehri/simulation/Mono/beta_chr1_GSE120610_Simulated.csv', data.table = F)
rownames(cd8s_1) <- cd8s_1$V1
cd8s_1$V1 <- NULL

monos_2 <- fread('/home/sepehri/simulation/Mono/beta_chr1_GSE131989.csv', data.table = F)
rownames(monos_2) <- monos_2$V1
monos_2$V1 <- NULL

cd8s_2 <- fread('/home/sepehri/simulation/Mono/beta_chr1_GSE131989_Simulated.csv', data.table = F)
rownames(cd8s_2) <- cd8s_2$V1
cd8s_2$V1 <- NULL

monos_3 <- fread('/home/sepehri/simulation/Mono/beta_chr1_GSE134429.csv', data.table = F)
rownames(monos_3) <- monos_3$V1
monos_3$V1 <- NULL

cd8s_3 <- fread('/home/sepehri/simulation/Mono/beta_chr1_GSE134429_Simulated.csv', data.table = F)
rownames(cd8s_3) <- cd8s_3$V1
cd8s_3$V1 <- NULL

monos_4 <- fread('/home/sepehri/simulation/Mono/beta_chr1_GSE184269.csv', data.table = F)
rownames(monos_4) <- monos_4$V1
monos_4$V1 <- NULL

cd8s_4 <- fread('/home/sepehri/simulation/Mono/beta_chr1_GSE184269_Simulated.csv', data.table = F)
rownames(cd8s_4) <- cd8s_4$V1
cd8s_4$V1 <- NULL

monos_0 <- fread("../data/monocytes/monocyte_raw_data/GSE56046_b.csv", data.table = F)
names(monos_0) <- monos_0[1,]
monos_0 <- monos_0[-1,]
rownames(monos_0) <- monos_0[,1]
monos_0 <- monos_0[,-1]

cd8s_0 <- fread('/home/sepehri/simulation/1202/S8_1202_Chr1.csv', data.table = F)
rownames(cd8s_0) <- cd8s_0$V1
cd8s_0$V1 <- NULL
names(cd8s_0) <- gsub("X", "", names(cd8s_0))

intersect_probes <- Reduce(intersect, list(rownames(monos_0), 
                                           rownames(monos_1), 
                                           rownames(monos_2), 
                                           rownames(monos_3), 
                                           rownames(monos_4)))

## filter on diff probes
diff_probes <- read.csv('../results/chr1_CD8-CD14_diff_probes.csv')
diff_probes$X <- NULL
intersect_probes <- intersect(intersect_probes, diff_probes$cluster_probes)

## select large data
large_sample_info <- read.csv('../data/monocytes/sample_info/large_age45to55.csv')
large_sample_info <- gsub("_.*", "", large_sample_info$title)
monos_0 <- monos_0[,large_sample_info]
cd8s_0 <- cd8s_0[,large_sample_info]

## check individual datasets
data <- fread("../data/monocytes/monocyte_raw_data/GSE56046_b.csv", data.table = F)
names(data) <- data[1,]
data <- data[-1,]
rownames(data) <- data[,1]
data <- data[,-1]

data <- data[intersect(diff_probes$cluster_probes, rownames(data)),]

test_monos <- read.csv('../data/monocytes/sample_info/test_age45to55.csv')
test_monos <- gsub("_.*", "", test_monos$title)
test_monos.B <- data[,test_monos]

names(test_monos.B) <- paste0("d0_mono_",1:dim(test_monos.B)[2])
sum(is.na(test_monos.B))
remove(data)

d1 <- fread("../data/monocytes/cd8/CD8_control_100.csv", data.table = F)
rownames(d1) <- d1$V1
d1$V1 <- NULL
names(d1) <- paste0("d1_CD8_",1:dim(d1)[2])
d1 <- d1[-unique(which(is.na(d1), arr.ind = T)[,1]),]
sum(is.na(d1))

d3 <- fread('~/../sepehri/simulation/final/28_CD8_CD14_control.csv', data.table = F)
names(d3)
rownames(d3) <- d3$V1
d3$V1<-NULL
names(d3) <- c(rep("d3_CD8", 28), rep("d3_mono", 28))
sum(is.na(d3))

intersect_probes <- Reduce(intersect, list(intersect_probes, 
                                           rownames(d3), 
                                           rownames(test_monos.B),
                                           rownames(d1)))
d3 <- d3[intersect_probes,]
d1 <- d1[intersect_probes,]
test_monos.B <- test_monos.B[intersect_probes,]

total <- cbind(d1, test_monos.B, d3)
sum(is.na(total))
total_tissues <- c(rep("CD8", dim(d1)[2]), rep("Monocyte", dim(test_monos.B)[2]), rep("CD8", 28), rep("Monocyte", 28))
total_tissues
total_batches <- c(rep("d1_CD8", dim(d1)[2]), rep("d0_mono", dim(test_monos.B)[2]),
                             rep("d3_CD8", 28), rep("d3_mono", 28))
total_batches
## select data
m <- monos_0
c <- cd8s_0

m <- m[intersect_probes,1:floor(dim(m)[2]/2)]
c <- c[intersect_probes,(floor(dim(c)[2]/2)+1):dim(c)[2]]

## merge data
all_m <- list(monos_0,monos_1,monos_2,monos_3,monos_4)
all_c <- list(cd8s_0, cd8s_1, cd8s_2, cd8s_3, cd8s_4)

m_tmp <- all_m[[1]]
c_tmp <- all_c[[1]]

m_tmp <- m_tmp[intersect_probes,1:floor(dim(m_tmp)[2]/2)]
c_tmp <- c_tmp[intersect_probes,(floor(dim(c_tmp)[2]/2)+1):dim(c_tmp)[2]]
print(dim(m_tmp))
print(dim(c_tmp))

m <- m_tmp
c <- c_tmp

for (i in 2:5){
  m_tmp <- all_m[[i]]
  c_tmp <- all_c[[i]]
  
  m_tmp <- m_tmp[intersect_probes,1:floor(dim(m_tmp)[2]/2)]
  c_tmp <- c_tmp[intersect_probes,(floor(dim(c_tmp)[2]/2)+1):dim(c_tmp)[2]]
  print(dim(m_tmp))
  print(dim(c_tmp))
  
  m <- cbind(m,m_tmp)
  c <- cbind(c,c_tmp)
}

## Plot individual Simulated Data
plot_df <- data.frame(t(cbind(m, c)))
tissue <- as.factor(c(rep("Monocyte", dim(m)[2]), rep("CD8", dim(c)[2])))
tissue
length(tissue)
batch <- as.factor(c(rep("0", floor(dim(monos_0)[2]/2)),
           rep("1", floor(dim(monos_1)[2]/2)),
           rep("2", floor(dim(monos_2)[2]/2)),
           rep("3", floor(dim(monos_3)[2]/2)),
           rep("4", floor(dim(monos_4)[2]/2)),
           rep("0", ceiling(dim(monos_0)[2]/2)),
           rep("1", ceiling(dim(monos_1)[2]/2)),
           rep("2", ceiling(dim(monos_2)[2]/2)),
           rep("3", ceiling(dim(monos_3)[2]/2)),
           rep("4", ceiling(dim(monos_4)[2]/2))))
length(batch)
pca_res <- prcomp(plot_df, scale. = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + theme_bw()
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point() + theme_bw()

## Plot Simulated data with Control data
tissue <- as.factor(total_tissues)
batch <- as.factor(total_batches)
mod <- model.matrix(~as.factor(tissue), data = total)
combat <- ComBat(dat=as.matrix(total), batch=as.factor(sub("_.*", "", batch)), mod=mod, par.prior=TRUE, prior.plots=FALSE)
dim(combat)

plot_df <- data.frame(t(combat))
pca_res <- prcomp(plot_df, scale. = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + ggtitle("CD8_Monocyte") + theme_bw()
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point() + ggtitle("CD8_Monocyte") + theme_bw()

selected_data <- cbind(combat, m, c) 
plot_df <- data.frame(t(selected_data))
tissue <- as.factor(c(total_tissues, rep("Monocyte", dim(m)[2]), rep("CD8", dim(c)[2])))
tissue
batch <- as.factor(c(total_batches, rep("sim_mono", dim(m)[2]), rep("sim_CD8", dim(c)[2])))
batch
pca_res <- prcomp(plot_df, scale. = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point() + theme_bw()
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + theme_bw()

## batch effect removal
total_data <- cbind(m, c)
tissue <- as.factor(c(rep("Monocyte", dim(m)[2]), rep("CD8", dim(c)[2])))
tissue
length(tissue)
batch <- as.factor(c(rep("0", floor(dim(monos_0)[2]/2)),
                     rep("1", floor(dim(monos_1)[2]/2)),
                     rep("2", floor(dim(monos_2)[2]/2)),
                     rep("3", floor(dim(monos_3)[2]/2)),
                     rep("4", floor(dim(monos_4)[2]/2)),
                     rep("0", ceiling(dim(monos_0)[2]/2)),
                     rep("1", ceiling(dim(monos_1)[2]/2)),
                     rep("2", ceiling(dim(monos_2)[2]/2)),
                     rep("3", ceiling(dim(monos_3)[2]/2)),
                     rep("4", ceiling(dim(monos_4)[2]/2))))
length(batch)
mod <- model.matrix(~as.factor(tissue), data = total_data)
combat <- ComBat(dat=as.matrix(total_data), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
dim(combat)
plot_df <- data.frame(t(combat))
pca_res <- prcomp(plot_df, scale. = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])
var_explained <- ((pca_res$sdev)^2/sum((pca_res$sdev)^2))*100
plot(var_explained)
n = 80
sum(var_explained[1:n])
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point() + theme_bw()
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + theme_bw()


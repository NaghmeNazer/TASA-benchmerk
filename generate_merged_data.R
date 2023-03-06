## Simulation datasets from different sources (GSE120610, GSE131989, GSE134429, GSE184269, GSE56046) are loaded
## to generate simulation scenario, different clusters of DMRs are chosen and simulation is done on them
## the generated dataset is gone through batch effect removal with combat from sva package
## simulation scenarios are saved for further analysis
library(data.table)
library(ggplot2)
library(sva)

simulate_cluster_probes <- function(m, c, p){
  p <- intersect(p, rownames(m))
  c_p <- m
  c_p[p,] <- c[p,]
  return(c_p)
}

SW_SB_cluster <- function(all_data, batch_col){
  batches = unique(batch_col)
  M = mean(rowMeans(all_data))
  SW = matrix(0, nrow = dim(all_data)[2], ncol = dim(all_data)[2])
  SB = SW
  for (b in batches){
    a <- all_data[batch_col == b, ]
    miu <- colMeans(a)
    a_miu <- as.matrix(sweep(a, 2, miu))
    SW <- SW + (t(a_miu) %*% a_miu)
    SB <- SB + dim(a)[1]*((miu - M) %*% t(miu-M))
  }
  tr.sb <- sum(diag(SB))
  tr.sw <- sum(diag(SW))
  dsc <- sqrt(tr.sb)/sqrt(tr.sw)
  return(list(tr.sb,tr.sw,dsc))
}

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
monos_4 <- monos_4[,-c(19,20)]

cd8s_4 <- fread('/home/sepehri/simulation/Mono/beta_chr1_GSE184269_Simulated.csv', data.table = F)
rownames(cd8s_4) <- cd8s_4$V1
cd8s_4$V1 <- NULL
cd8s_4 <- cd8s_4[,-c(19,20)]

monos_0 <- fread("../data/monocytes/monocyte_raw_data/GSE56046_b.csv", data.table = F)
names(monos_0) <- monos_0[1,]
monos_0 <- monos_0[-1,]
rownames(monos_0) <- monos_0[,1]
monos_0 <- monos_0[,-1]

cd8s_0 <- fread('/home/sepehri/simulation/1202/S8_1202_Chr1.csv', data.table = F)
rownames(cd8s_0) <- cd8s_0$V1
cd8s_0$V1 <- NULL
names(cd8s_0) <- gsub("X", "", names(cd8s_0))

## select large data
large_sample_info <- read.csv('data/sample_info/large_age45to55.csv')
large_sample_info <- gsub("_.*", "", large_sample_info$title)
monos_0 <- monos_0[,large_sample_info]
cd8s_0 <- cd8s_0[,large_sample_info]

## intersect probes
intersect_probes <- Reduce(intersect, list(rownames(monos_0), 
                                           rownames(monos_1), 
                                           rownames(monos_2), 
                                           rownames(monos_3), 
                                           rownames(monos_4)))
bad_qual_probes <- fread('/home/sepehri/simulation/Mono/Unwanted_Probes_chr1.csv', data.table = F)
bad_qual_probes <- bad_qual_probes$V1
intersect_probes <- setdiff(intersect_probes, bad_qual_probes)
sum(bad_qual_probes %in% intersect_probes)

monos_0 <- monos_0[intersect_probes,1:14]
sum(is.na(monos_0))
cd8s_0 <- cd8s_0[intersect_probes,1:14]
sum(is.na(cd8s_0))

monos_1 <- monos_1[intersect_probes,1:6]
sum(is.na(monos_1))
cd8s_1 <- cd8s_1[intersect_probes,1:6]
sum(is.na(cd8s_1))

monos_2 <- monos_2[intersect_probes,1:4]
sum(is.na(monos_2))
cd8s_2 <- cd8s_2[intersect_probes,1:4]
sum(is.na(cd8s_2))

monos_3 <- monos_3[intersect_probes,1:4]
sum(is.na(monos_3))
cd8s_3 <- cd8s_3[intersect_probes,1:4]
sum(is.na(cd8s_3))

monos_4 <- monos_4[intersect_probes,1:4]
sum(is.na(monos_4))
cd8s_4 <- cd8s_4[intersect_probes,1:4]
sum(is.na(cd8s_4))

## read cluster probes
cluster_probes <- read.csv('../results/chr1_cluster_probes.csv')
cluster_probes$X <- NULL
cluster_probes <- cluster_probes[cluster_probes$cluster_probes %in% intersect_probes,]
dim(cluster_probes)[1]

## generate merged data -- no simulation
all_m <- list(monos_0,monos_1,monos_2,monos_3,monos_4)
all_c <- list(monos_0, monos_1, monos_2, monos_3, monos_4)

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

## Plot merged data
total_data <- cbind(m, c)
plot_df <- data.frame(t(total_data))
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

SW_SB_cluster(total_data, tissue)
SW_SB_cluster(total_data, batch)

## save data
fwrite(total_data, "../data/small_simulated/small_data_raw_no_simulation_batch_not_corrected.csv", col.names = T, row.names = T)


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
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point() + theme_bw()
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + theme_bw()
SW_SB_cluster(combat, tissue)
SW_SB_cluster(combat, batch)

## save data
combat <- data.frame(combat)
rownames(combat)
names(combat) <- names(total_data)
fwrite(combat, "../data/small_simulated/small_data_raw_no_simulation_batch_corrected.csv", col.names = T, row.names = T)


## simulate cluster probes
cd8s_p_0 <- simulate_cluster_probes(monos_0, cd8s_0, cluster_probes$cluster_probes)
cd8s_p_1 <- simulate_cluster_probes(monos_1, cd8s_1, cluster_probes$cluster_probes)
cd8s_p_2 <- simulate_cluster_probes(monos_2, cd8s_2, cluster_probes$cluster_probes)
cd8s_p_3 <- simulate_cluster_probes(monos_3, cd8s_3, cluster_probes$cluster_probes)
cd8s_p_4 <- simulate_cluster_probes(monos_4, cd8s_4, cluster_probes$cluster_probes)

## generate merged data
all_m <- list(monos_0,monos_1,monos_2,monos_3,monos_4)
all_c <- list(cd8s_p_0, cd8s_p_1, cd8s_p_2, cd8s_p_3, cd8s_p_4)

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

## Plot merged data
total_data <- cbind(m, c)
plot_df <- data.frame(t(total_data))
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
total_data_info <- data.frame("Sample_Name" = names(total_data), "Tissue" = tissue, "Batch" = batch)
pca_res <- prcomp(plot_df, scale. = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + theme_bw()
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point() + theme_bw()

SW_SB_cluster(total_data, tissue)
SW_SB_cluster(total_data, batch)

## save data
fwrite(total_data, "../data/small_simulated/0.4_small_data_raw_batch_not_corrected.csv", col.names = T, row.names = T)

fwrite(total_data_info, "../data/small_simulated/small_data_info.csv", col.names = T, row.names = T)
b <- fread("data/sample_info/small_data_info.csv", data.table = F)
b <- b[,-1]

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
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point() + theme_bw()
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + theme_bw()
SW_SB_cluster(combat, tissue)
SW_SB_cluster(combat, batch)

## save data
combat <- data.frame(combat)
rownames(combat)
names(combat) <- names(total_data)
fwrite(combat, "../data/small_simulated/0.4_small_data_raw_batch_corrected.csv", col.names = T, row.names = T)





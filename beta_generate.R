# Preprocess and generate beta value for large size data
# Make three separate matrices of Methylated signal intensities, Unmethylated signal intensities, and P-values
# You will need rs_af0.05.csv file (probes with mutation more than 5% allele freq),
# and non_specific_sites.csv file (probes attaching to multiple loci in genome)
library(data.table)
library(ggvenn)

intensities.to.beta <- function(M, U, alpha = 100){
  B <- M/(M+U+alpha)
  return(B)
}

bad.detected.probes <- function(P){
  return(unique(rownames(which(P > 0.05, arr.ind = T))))
}

negative.probes <- function(M, U){
  return(union(unique(rownames(which(M < 0, arr.ind = T))), unique(rownames(which(U < 0, arr.ind = T)))))
}

rs_probes <- function(){
  rs <- read.csv('data/monocytes/manifest/rs_af0.05.csv')
  return(unique(rs$X))
}

non_specific.probes <- function(){
  ps <- read.csv('/media/mehrmohammadi_hdd/nazer/data/monocytes/manifest/non_specific_sites.csv')
  return(unique(ps[,2]))
}

deleted.probes <- function(M,U,P){
  high.p.probes <- bad.detected.probes(P)
  neg.probes <- negative.probes(M, U)
  nspec.probes <- non_specific.probes()
  snp.probes <- rs_probes()
  deleted.probes <- unique(c(high.p.probes, neg.probes, nspec.probes, snp.probes))
  return(deleted.probes)
}

preprocess_to_beta <- function(M, U, P){
  bad.probes <- deleted.probes(M,U,P)
  B <- intensities.to.beta(M,U)
  B <- B[!(rownames(B) %in% bad.probes),]
  return(B)
}

P <- fread('../data/monocytes/monocyte_raw_data/GSE56046_p.txt', data.table = F)
rownames(P) <- P$ID_REF
P$ID_REF <- NULL
names(P) <- gsub(".detectionPval", "", names(P))
M <- fread('../data/monocytes/monocyte_raw_data/GSE56046_m.txt', data.table = F)
rownames(M) <- rownames(P)
M$ID_REF <- NULL
names(M) <- names(P)
U <- fread('../data/monocytes/monocyte_raw_data/GSE56046_u.txt', data.table = F)
rownames(U) <- rownames(P)
U$ID_REF <- NULL
names(U) <- names(P)

#bad.probes <- deleted.probes(M,U,P)
#names(bad.probes) <- c("High-P", "Negative", "Non-specific", "SNP")

#ggvenn( bad.probes, 
#  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#  stroke_size = 0.5, set_name_size = 4)

N <- 50000
P1 <- P[1:N,]
P1 <- P1[complete.cases(P1),]
M1 <- M[1:N,]
M1 <- M1[complete.cases(M1),]
U1 <- U[1:N,]
U1 <- U1[complete.cases(U1),]
B <- preprocess_to_beta(M1, U1, P1)
remove(M1)
remove(U1)
remove(P1)
P <- P[-c(1:N),]
M <- M[-c(1:N),]
U <- U[-c(1:N),]
gc()

while (dim(P)[1] > 0){
  P1 <- P[1:N,]
  P1 <- P1[complete.cases(P1),]
  M1 <- M[1:N,]
  M1 <- M1[complete.cases(M1),]
  U1 <- U[1:N,]
  U1 <- U1[complete.cases(U1),]
  B1 <- preprocess_to_beta(M1, U1, P1)
  B <- rbind(B, B1)
  remove(M1)
  remove(U1)
  remove(P1)
  remove(B1)
  P <- P[-c(1:N),]
  M <- M[-c(1:N),]
  U <- U[-c(1:N),]
  print("on-track")
  gc()
}

remove(M)
remove(P)
remove(U)
gc()

B <- round(B, digits = 2)
fwrite(B, "../data/monocytes/monocyte_raw_data/GSE56046_b_rounded.csv", col.names = T, row.names = T)


## read beta values
b <- fread("../data/monocytes/monocyte_raw_data/GSE56046_b.csv", data.table = F)
names(b) <- b[1,]
b <- b[-1,]
rownames(b) <- b[,1]
b <- b[,-1]
bb <- b[1:10,1:10]
bad.probes <- unique(c(non_specific.probes(), rs_probes()))


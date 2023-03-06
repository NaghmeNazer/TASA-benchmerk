## The SVM classifier to check the simulation datasets
## simulated data set alongside real samples of monocytes, CD8, and breast tissue are loaded
## the real data is splitted to test and train
## first 10 PCS of trained data is used and SVM is trained with 10-fold cross validation
## test and simulated data is projected to PCA space and classified with the trained model
## PCA plots are generated
library(data.table)
library(caTools)
library(e1071)
library(ggplot2)
library(varhandle)
library(matrixStats)
library(limma)
library(sva)
library(pROC)
library(caret)

#### functions
# make beta value
intensities.to.beta <- function(M, U, alpha = 100){
  B <- M/(M+U+alpha)
  return(B)
}

# find bad detected probes
bad.detected.probes <- function(P){
  return(unique(rownames(which(P > 0.05, arr.ind = T))))
}

# fine probes with negative intensity
negative.probes <- function(M, U){
  return(union(unique(rownames(which(M < 0, arr.ind = T))), unique(rownames(which(U < 0, arr.ind = T)))))
}

# preprocess (probe qc) and generate beta values
preprocess_to_beta <- function(M, U, P){
  high.p.probes <- bad.detected.probes(P)
  neg.probes <- negative.probes(M, U)
  B <- intensities.to.beta(M,U)
  B <- B[!(rownames(B) %in% c(high.p.probes, neg.probes)),]
  return(B)
}

### read datasets

#data <- fread("../data/monocytes/monocyte_raw_data/GSE56046_b.csv", data.table = F)
data <- fread("/media/mehrmohammadi_hdd/nazer/data/monocytes/monocyte_raw_data/GSE56046_b.csv", data.table = F)
names(data) <- data[1,]
data <- data[-1,]
rownames(data) <- data[,1]
data <- data[,-1]

test_monos <- read.csv('/media/mehrmohammadi_hdd/nazer/data/monocytes/sample_info/test_age45to55.csv')
test_monos <- gsub("_.*", "", test_monos$title)
test_monos.B <- data[,test_monos]

names(test_monos.B) <- paste0("d0_mono_",1:dim(test_monos.B)[2])
sum(is.na(test_monos.B))
remove(data)

d1 <- fread("/media/mehrmohammadi_hdd/nazer/data/monocytes/cd8/CD8_control_100.csv", data.table = F)
rownames(d1) <- d1$V1
d1$V1 <- NULL
names(d1) <- paste0("d1_CD8_",1:dim(d1)[2])
d1 <- d1[-unique(which(is.na(d1), arr.ind = T)[,1]),]
sum(is.na(d1))

d3 <- fread('/home/psepehri/28_CD8_CD14_control.csv', data.table = F)
names(d3)
rownames(d3) <- d3$V1
d3$V1<-NULL
names(d3) <- c(rep("d3_CD8", 28), rep("d3_mono", 28))
sum(is.na(d3))

d4 <- fread('/media/mehrmohammadi_hdd/sepehri/sepehri/sepehri/simulation/Mono_Breast/beta_breast_control.csv', data.table = F)
rownames(d4) <- d4$V1
d4$V1 <- NULL
names(d4) <- paste0("d4_Breast_",1:dim(d4)[2])
sum(is.na(d4))

## intersect probes of all datasets
intersect_probes <- Reduce(intersect, list(rownames(test_monos.B),
                                           rownames(d1),
                                           rownames(d3),
                                           rownames(d4)))
d3 <- d3[intersect_probes,]
d1 <- d1[intersect_probes,]
test_monos.B <- test_monos.B[intersect_probes,]
d4 <- d4[intersect_probes,]

total <- cbind(d1, test_monos.B, d3, d4)
sum(is.na(total))

CD8_simulated_data <- fread('/home/psepehri/S10_1202_Chr1.csv', header = T, data.table = F)
rownames(CD8_simulated_data) <- CD8_simulated_data$V1
CD8_simulated_data$V1 <- NULL
names(CD8_simulated_data) <- gsub("X", "", names(CD8_simulated_data))
CD8_simulated_data <- CD8_simulated_data[,!(names(CD8_simulated_data) %in% test_monos)]
names(CD8_simulated_data) <- paste0("sim_CD8_", 1:1092)

Breast_simulated_data <- fread('/home/psepehri/beta_chr1_1202_breast.csv', header = T, data.table = F)
rownames(Breast_simulated_data) <- Breast_simulated_data$V1
Breast_simulated_data$V1 <- NULL
names(Breast_simulated_data) <- gsub("X", "", names(Breast_simulated_data))
Breast_simulated_data <- Breast_simulated_data[,!(names(Breast_simulated_data) %in% test_monos)]
names(Breast_simulated_data) <- paste0("sim_Breast_", 1:1092)

final_probes <- Reduce(intersect, list(intersect_probes,
                                       rownames(CD8_simulated_data),
                                       rownames(Breast_simulated_data)))

CD8_simulated_data <- CD8_simulated_data[final_probes,]
Breast_simulated_data <- Breast_simulated_data[final_probes,]
simulated_data <- cbind(CD8_simulated_data, Breast_simulated_data)
tail(names(simulated_data))
sum(is.na(simulated_data))

total <- total[final_probes,]
total_sim <- cbind(total, simulated_data)
names(total_sim)

## PCA plot of datasets
plot_df <- total_sim
sum(is.na(total_sim))
tissue <- rep("CD8", dim(total_sim)[2])
tissue[grepl("mono", names(total_sim))] <- "mono"
tissue[grepl("Breast", names(total_sim))] <- "Breast"
tissue <- as.factor(tissue)
tissue
batch <- as.factor(c(rep("d1_cd8", dim(d1)[2]), rep("d0_mono", dim(test_monos.B)[2]), names(d3), rep("d4_breast", dim(d4)[2]), rep("sim_cd8", 1092), rep("sim_breast", 1092)))
batch
pca_res <- prcomp(data.frame(t(plot_df)), scale. = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])

pdf("../results/plots/tri_total_sim.pdf")
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + ggtitle("CD8_Monocyte_Breast") + theme_bw()
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point() + ggtitle("CD8_Monocyte_Breast") + theme_bw()
dev.off()

plot_df <- total
tissue <- rep("CD8", dim(total)[2])
tissue[grepl("mono", names(total))] <- "mono"
tissue[grepl("Breast", names(total))] <- "Breast"
tissue <- as.factor(tissue)
tissue
batch <- as.factor(c(rep("d1_cd8", dim(d1)[2]), rep("d0_mono", dim(test_monos.B)[2]), names(d3),  rep("d4_breast", dim(d4)[2])))
batch
pca_res <- prcomp(data.frame(t(plot_df)), scale. = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])

pdf("../results/plots/tri_total.pdf")
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + ggtitle("CD8_Monocyte_Breast") + theme_bw()
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point() + ggtitle("CD8_Monocyte_Breast") + theme_bw()
dev.off()
#plot_df <- total_sim[,101:414]
#tissue <- tissue[101:414]
#batch <- batch[101:414]
#pca_res <- prcomp(data.frame(t(plot_df)), scale. = TRUE)
#pca_res.x <- data.frame(pca_res[["x"]][,1:2])
#ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + ggtitle("CD8_Monocyte") + theme_bw()


## batch effect removal
tissue <- rep("CD8", dim(total)[2])
tissue[grepl("mono", names(total))] <- "mono"
tissue[grepl("Breast", names(total))] <- "Breast"
tissue <- as.factor(tissue)
tissue
batch <- as.factor(c(rep("d1_cd8", dim(d1)[2]), rep("d0_mono", dim(test_monos.B)[2]), names(d3),  rep("d4_breast", dim(d4)[2])))
batch
mod <- model.matrix(~as.factor(tissue), data = total)
combat <- ComBat(dat=as.matrix(total), batch=as.factor(sub("_.*", "", batch)), mod=mod, par.prior=TRUE, prior.plots=FALSE)
dim(combat)
plot_df <- data.frame(t(combat))
pca_res <- prcomp(plot_df, scale. = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + ggtitle("CD8_Monocyte") + theme_bw()

## PCA plot of batch effect removed data
combat_sim <- cbind(combat, simulated_data)
names(combat_sim)
plot_df <- combat_sim
sum(is.na(combat_sim))
tissue <- rep("CD8", dim(combat_sim)[2])
tissue[grepl("mono", names(combat_sim))] <- "mono"
tissue <- as.factor(tissue)
tissue
batch <- as.factor(c(rep("d1_cd8", dim(d1)[2]), rep("d0_mono", dim(test_monos.B)[2]), names(d3), rep("sim_cd8", 1092)))
batch
pca_res <- prcomp(data.frame(t(plot_df)), scale. = TRUE)
pca_res.x <- data.frame(pca_res[["x"]][,1:2])
ggplot(pca_res.x, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point() + ggtitle("CD8_Monocyte") + theme_bw()

## Training
df <- data.frame(t(total))
#df <- data.frame(t(combat))
df$y <- "CD8"
df$y[grepl("mono", rownames(df))] <- "mono"
df$y[grepl("Breast", rownames(df))] <- "Breast"
df$y <- as.factor(df$y)
df$y

## split test and train samples
set.seed(123)
split <- sample.split(df$y, SplitRatio = 0.8)
training_set = subset(df, split == TRUE)
test_set = subset(df, split == FALSE)

## generate PCA transform for training samples
pca_res <- prcomp(training_set[-dim(training_set)[2]], scale. = TRUE)
var_explained <- ((pca_res$sdev)^2/sum((pca_res$sdev)^2))*100
#plot(var_explained)
n = 10
sum(var_explained[1:n])

training_set_pcs <- scale(training_set[-dim(training_set)[2]], pca_res$center, pca_res$scale) %*% pca_res$rotation 
training_set_pcs <- data.frame(training_set_pcs[,1:n])
training_set_pcs$y <- training_set$y
training_set_pcs[-(n+1)] = scale(training_set_pcs[-(n+1)])

#training_set_pcs <- training_set[,c(boruta_signif, "y")]

## learn SVM classifier with 10-fold cv
classifier = svm(formula = y ~ .,
                 data = training_set_pcs,
                 cross=10,
                 type = 'C-classification',
                 kernel = 'linear')

## get predictions of training samples and compute decision values
y_train_pred = predict(classifier, newdata = training_set_pcs[-dim(training_set_pcs)[2]], decision.values = T)
cm_training = table(training_set_pcs$y, y_train_pred)
cm_training
dec.val_training <- mean(abs(attr(y_train_pred, "decision.values")))
dec.val_training

## get predictions of test samples and compute decision values
test_set_pcs <- scale(test_set[-dim(test_set)[2]], pca_res$center, pca_res$scale) %*% pca_res$rotation 
test_set_pcs <- data.frame(test_set_pcs[,1:n])
test_set_pcs$y <- test_set$y
test_set_pcs[-(n+1)] = scale(test_set_pcs[-(n+1)])

y_test_pred = predict(classifier, newdata = test_set_pcs[-dim(test_set_pcs)[2]], decision.values = T)
cm_test = table(test_set_pcs$y, y_test_pred)
cm_test
dec.val_test <- mean(abs(attr(y_test_pred, "decision.values")))
dec.val_test

## test simulated data
## read simulation data (monocyte, CD8, breast)
data <- fread("/media/mehrmohammadi_hdd/nazer/data/monocytes/monocyte_raw_data/GSE56046_b.csv", data.table = F)
names(data) <- data[1,]
data <- data[-1,]
rownames(data) <- data[,1]
data <- data[,-1]
data <- data[rownames(total),]
#data <- data[rownames(combat),]

CD8_simulated_data <- fread('/home/psepehri/S10_1202_Chr1.csv', header = T, data.table = F)
rownames(CD8_simulated_data) <- CD8_simulated_data$V1
CD8_simulated_data$V1 <- NULL
names(CD8_simulated_data) <- gsub("X", "", names(CD8_simulated_data))
CD8_simulated_data <- CD8_simulated_data[,!(names(CD8_simulated_data) %in% test_monos)]
CD8_simulated_data <- CD8_simulated_data[rownames(total),]
#CD8_simulated_data <- CD8_simulated_data[rownames(combat),]
CD8_simulated_data <- cbind(data[names(CD8_simulated_data)[1:546]], CD8_simulated_data[547:1092])
sum(is.na(CD8_simulated_data))
names(CD8_simulated_data) <- c(paste0("sim_Monocyte_", 1:546), paste0("sim_CD8_", 547:1092))

Breast_simulated_data <- fread('/home/psepehri/beta_chr1_1202_breast.csv', header = T, data.table = F)
rownames(Breast_simulated_data) <- Breast_simulated_data$V1
Breast_simulated_data$V1 <- NULL
names(Breast_simulated_data) <- gsub("X", "", names(Breast_simulated_data))
Breast_simulated_data <- Breast_simulated_data[,!(names(Breast_simulated_data) %in% test_monos)]
Breast_simulated_data <- Breast_simulated_data[rownames(total),]
#CD8_simulated_data <- CD8_simulated_data[rownames(combat),]
Breast_simulated_data <- cbind(data[names(Breast_simulated_data)[1:546]], Breast_simulated_data[547:1092])
sum(is.na(Breast_simulated_data))
names(Breast_simulated_data) <- c(paste0("sim_Monocyte_", 1:546), paste0("sim_Breast_", 547:1092))

simulated_data <- cbind(CD8_simulated_data, Breast_simulated_data[,547:1092])
simulated_data <- data.frame(t(simulated_data))
rownames(simulated_data)
simulated_data$y <- c(rep("Monocyte", 546), rep("CD8", 546), rep("Breast", 546))
simulated_data$y <- as.factor(simulated_data$y)
simulated_data$y

## project simulated data on PCA space
simulated_data_pcs <- scale(simulated_data[-dim(simulated_data)[2]], pca_res$center, pca_res$scale) %*% pca_res$rotation
#N = 204
#simulated_data_pcs <- scale(simulated_data[103:204,-dim(simulated_data)[2]], pca_res$center, pca_res$scale) %*% pca_res$rotation 
dim(simulated_data_pcs)
simulated_data_pcs <- data.frame(simulated_data_pcs[,1:n])
simulated_data_pcs$y <- simulated_data$y
simulated_data_pcs[-(n+1)] = scale(simulated_data_pcs[-(n+1)])

## plot simulated data and real samples (test and train datasets)
plot_df <- rbind(training_set_pcs[,c(1,2,n+1)], test_set_pcs[,c(1,2,n+1)], simulated_data_pcs[,c(1,2,n+1)])
plot_df$y[plot_df$y == "mono"] <- "Monocyte"
tissue <- as.factor(as.character(plot_df$y))
batch <- as.factor(c(rep("Training", dim(training_set_pcs)[1]), 
                     rep("Test", dim(test_set_pcs)[1]),
                     rep("Simulated", dim(simulated_data_pcs)[1])))
pdf("../nazer/results/plots/PCA-SVM-plt.pdf")
ggplot(plot_df, aes(x=PC1, y=PC2, colour = batch, shape = tissue)) + geom_point()  + theme_bw() + theme(text = element_text(size = 20))
#ggplot(plot_df, aes(x=PC1, y=PC2, colour = tissue, shape = batch)) + geom_point()  + theme_bw()
dev.off()

## get predictions of simulated samples and compute decision values
y_sim_pred = predict(classifier, newdata = simulated_data_pcs, decision.values = T)
cm_sim = table(simulated_data_pcs$y, y_sim_pred)
cm_sim
errors = c(which(y_sim_pred[1:546] == "CD8"), which(y_sim_pred[547:1092] == "mono"))
y_sim_pred.decval <- attr(y_sim_pred, "decision.values")
dec.val_sim <- (sum(abs(y_sim_pred.decval)) - 2*sum(abs(y_sim_pred.decval[errors])))/length(y_sim_pred.decval)
dec.val_sim

## ROC analysis
roc_svm_train <- roc(response = training_set_pcs$y, predictor =as.numeric(y_train_pred))
roc_svm_test <- roc(response = test_set_pcs$y, predictor =as.numeric(y_test_pred))
levels(simulated_data_pcs$y) <- c("CD8", "mono", "Breast")
roc_svm_simulated <- roc(response = simulated_data_pcs$y, predictor =as.numeric(y_sim_pred))

best.point <- coords(roc_svm_train, "best", ret=c("sensitivity", "specificity"), transpose = FALSE)

roc.legend <- paste0(sprintf('Trainig Set AUC = %.2f', roc_svm_train$auc), "\n",
                     sprintf('Test Set AUC = %.2f', roc_svm_test$auc), "\n",
                     sprintf('Simulated Set AUC = %.2f', roc_svm_simulated$auc))
ggroc(list(Training= roc_svm_train, Test= roc_svm_test, Simulated = roc_svm_simulated), legacy.axes = TRUE) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed") + 
  annotate("text",  x=0.6, y = 0, hjust = 'left', vjust = 'bottom', label = roc.legend, vjust=1.25, hjust=-0.125, size = 3) +
  ggtitle("ROC Title") + theme_bw()


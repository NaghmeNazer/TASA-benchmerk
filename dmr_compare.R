### Script to compare different dmr finding methods applied to datasets
setwd("code/")
library(data.table)
library(dplyr)

###############
### Functions
### function to map probes to their genomic loci
probes_to_location <- function(probes){
  manifest <- fread("/media/mehrmohammadi_hdd/nazer/data/monocytes/manifest/450k_manifest.csv", data.table = F)
  rownames(manifest) <- manifest$IlmnID
  probes$CHR <- manifest[rownames(probes), "CHR"]
  probes$MAPINFO <- manifest[rownames(probes), "MAPINFO"]
  probes_gpd <- probes %>% group_by(cluster_numbers)
  probes_loc <- data.frame(probes_gpd %>% dplyr::summarise(chr = unique(CHR), start = min(MAPINFO), end = max(MAPINFO)+1, length = end-start+1, nprobe = length(cluster_probes)))
  probes_loc <- probes_loc[order(probes_loc$start),]
  return(probes_loc)
}

### function that generates a binary signal of DMR patterns.
### loci inside a dmr will be 1, other loci will be 0
make_signal <- function(dmrs, cluster_number = F){
  if (dim(dmrs)[1] == 0){
    signal <- rep(0,1)
  }
  else{
    max <- max(dmrs$end)
    signal <- rep(0, max)
    for (i in 1:dim(dmrs)[1]){
      if(cluster_number){
        signal[dmrs[i, "start"]:dmrs[i,"end"]] <- dmrs[i, "cluster_numbers"]
      }
      else{
        signal[dmrs[i, "start"]:dmrs[i,"end"]] <- 1
      }
    }
  }
  return(signal)
}

### function to compute intersection of binary signals
intersect_signal <- function(s1, s2){
  max <- max(length(s1), length(s2))
  if (length(s1) < max){
    s1 <- c(s1, rep(0, max - length(s1)))
  }
  if (length(s2) < max){
    s2 <- c(s2, rep(0, max - length(s2)))
  }
  s_i <- s1*s2
  return(s_i)
}

### function to calculate the overlap cverage in DMRs
compute_coverage <- function(dmr, s_i){
  if(dim(dmr)[1] == 0){
    return(dmr)
  }
  else{
    dmr$covered <- NA
    for (i in 1:dim(dmr)[1]){
      dmr[i, "covered"] <- sum(s_i[dmr[i, "start"]:dmr[i, "end"]])
    }
    dmr$length <- dmr$end - dmr$start + 1
    dmr$covered_ratio <- (dmr$covered)/(dmr$length)
    return(dmr)
  }
}

### function that adds evaluation metrics to the input dataframe
### Evaluation metrics are TP, FP, FN, Recall, Precision
evaluation_metrics <- function(f, evaluation, true_dmr, result_dmr){
  evaluation[f, "File"] <- f
  evaluation[f, "TP"] <- sum(true_dmr$TP)
  
  evaluation[f, "FP"] <- sum(result_dmr$FP)
  evaluation[f, "FN"] <- sum(true_dmr$FN)
  evaluation[f, "Recall"] <- evaluation[f, "TP"]/(evaluation[f, "TP"]+evaluation[f, "FN"])
  evaluation[f, "Precision"] <- evaluation[f, "TP"]/(evaluation[f, "TP"]+evaluation[f, "FP"])
  return(evaluation)
}

################
### Main
### list outputs of different DMR tools stored in a directory
dir = "../results/combp_DMRs/small_CD8/"
files <- list.files(dir)
files <- files[!(grepl("betaQN", files) | grepl("evaluation", files))]
files

### make avaluation dataframe and fill each column as needed
evaluation <- data.frame(matrix(nrow = length(files), ncol = 8))
names(evaluation) <- c("Scenario", "File", "Method", "TP", "FP", "FN", "Recall", "Precision")
rownames(evaluation) <- files
#evaluation$Method <- gsub(".csv","",gsub("^.*\\_","",rownames(evaluation)))
evaluation$Method <- substr(dir,12,16)
evaluation$Scenario <- gsub("^.*\\/","",substr(dir,1,nchar(dir)-1))

### read average of probe differences in DMRs when you want to compare tools in high/low differentiated regins
if(grepl("CD8", evaluation$Scenario[1])){
  probe_diff <- read.csv('/media/mehrmohammadi_hdd/sepehri/sepehri/sepehri/simulation/CD8_CD14_difference.csv')
  names(probe_diff) <- c('cluster_probes', 'diff')
}
if(grepl("Breast", evaluation$Scenario[1])){
  probe_diff <- read.csv('/media/mehrmohammadi_hdd/sepehri/sepehri/sepehri/simulation/Mono_Breast/breast_CD14_difference.csv')
  probe_diff$chr <- NULL
  names(probe_diff) <- c('cluster_probes', 'diff')
}


coverage_threshod <- 0.2
intensity_threshold <- 0.3

evaluation_high <- evaluation
evaluation_low <- evaluation

### compute evaluation metrics with predefined functions
### the metrics are calculated for all DMR tools
### differences in the format of outputs are handled in this for loop
for (f in files){
  print(f)
  result_dmr <- fread(paste0(dir,f), data.table = F)
  #result_dmr <- fread(f, data.table = F)
  names(result_dmr) <- gsub("^.*\\.","",names(result_dmr))
  
  names(result_dmr)[names(result_dmr) == 'stop'] <- "end"
  if(evaluation[f, "Method"] == "DMRcate"){
    result_dmr$coord <- gsub(".*:", "", result_dmr$coord)
    result_dmr$end <- as.numeric(gsub(".*-", "", result_dmr$coord))
    result_dmr$start <- as.numeric(gsub("-.*", "", result_dmr$coord))
  }
  result_dmr <- result_dmr[order(result_dmr$start),]
  result_dmr$length <- NULL
  result_dmr$length <- result_dmr$end - result_dmr$start + 1
  
  names(result_dmr)[names(result_dmr) == 'HMFDR'] <- "fdr"
  names(result_dmr)[names(result_dmr) == 'fwer'] <- "fdr"
  if (!(evaluation[f, "Method"] == "ProbeLasso" | evaluation[f, "Method"] == "Cohcap")){
    result_dmr <- result_dmr[result_dmr$fdr <= 0.05,]
  }
  result_dmr_signal <- make_signal(result_dmr)
  
  if (grepl("0.1", f)){
    probes <- fread('../results/chr1_0.1_cluster_probes.csv', data.table = F)
    rownames(probes) <- probes$cluster_probes
    probes$V1 <- NULL
    true_dmr <- probes_to_location(probes)
    true_dmr_signal <- make_signal(true_dmr)
    true_dmr_signal_cs <- make_signal(true_dmr, T)
  }
  if (grepl("0.2", f)){
    probes <- fread('../results/chr1_0.2_cluster_probes.csv', data.table = F)
    rownames(probes) <- probes$cluster_probes
    probes$V1 <- NULL
    true_dmr <- probes_to_location(probes)
    true_dmr_signal <- make_signal(true_dmr)
    true_dmr_signal_cs <- make_signal(true_dmr, T)
  }
  if (grepl("0.4", f)){
    probes <- fread('../results/chr1_cluster_probes.csv', data.table = F)
    rownames(probes) <- probes$cluster_probes
    probes$V1 <- NULL
    true_dmr <- probes_to_location(probes)
    true_dmr_signal <- make_signal(true_dmr)
    true_dmr_signal_cs <- make_signal(true_dmr, T)
  }
  if (grepl("no_simulation", f) | grepl("ns", f)){
    true_dmr <- data.frame(matrix(nrow = 0, ncol = 6))
    names(true_dmr) <- c("cluster_numbers", "chr", "start", "end", "length", "nprobe")
    true_dmr_signal <- make_signal(true_dmr)
    true_dmr_signal_cs <- make_signal(true_dmr, T)
  }
  
  if (!(grepl("no_simulation", f) | grepl("ns", f))){
    cluster_diff <- merge(probe_diff, probes, by = "cluster_probes")
    clusters_diff_gpd <- cluster_diff %>% group_by(cluster_numbers)
    clusters_diff_mean <- data.frame(clusters_diff_gpd %>% dplyr::summarise(mean_diff = mean(diff)))
  }
  i_signal <- intersect_signal(result_dmr_signal, true_dmr_signal)
  i_signal_cs <- intersect_signal(result_dmr_signal, true_dmr_signal_cs)
  
  true_dmr <- compute_coverage(true_dmr, i_signal)
  true_dmr$FN <- true_dmr$covered_ratio < coverage_threshod
  true_dmr$TP <- true_dmr$covered_ratio >= coverage_threshod
  if (!(grepl("no_simulation", f) | grepl("ns", f))){
    true_dmr <- merge(true_dmr, clusters_diff_mean, by = "cluster_numbers")
  }
  
  result_dmr <- compute_coverage(result_dmr, i_signal)
  #result_dmr$TP <- result_dmr$covered_ratio >= coverage_threshod
  result_dmr$FP <- result_dmr$covered_ratio < coverage_threshod
  
  if (!(grepl("no_simulation", f) | grepl("ns", f))){
    result_dmr$mean_diff <- NA
    
    for (r in 1:dim(result_dmr)[1]){
      cs <- unique(i_signal_cs[result_dmr[r, "start"]:result_dmr[r,"end"]])
      cs <- cs[cs!=0]
      if (length(cs) == 0){
        result_dmr[r, "mean_diff"] <- 0
      }
      else{
        for (c in cs){
          result_dmr[r, "mean_diff"] <- mean(clusters_diff_mean[clusters_diff_mean$cluster_numbers %in% cs, "mean_diff"])
        }
      }
    }
  }
  
  evaluation <- evaluation_metrics(f, evaluation, true_dmr, result_dmr)
  if (!(grepl("no_simulation", f) | grepl("ns", f))){
    evaluation_high <- evaluation_metrics(f, evaluation_high, true_dmr[abs(true_dmr$mean_diff) >= intensity_threshold,], result_dmr[abs(result_dmr$mean_diff) >= intensity_threshold,])
    evaluation_low <- evaluation_metrics(f, evaluation_low, true_dmr[abs(true_dmr$mean_diff) < intensity_threshold,], result_dmr[abs(result_dmr$mean_diff) < intensity_threshold,])
  }
  
  rm(cluster_diff, clusters_diff_gpd, clusters_diff_mean, probes,
     true_dmr, true_dmr_signal, true_dmr_signal_cs, i_signal, i_signal_cs, result_dmr, result_dmr_signal)
  gc()
}


### Write results
#fwrite(evaluation, paste0("../results/sepehri_DMRs/",gsub("^.*\\/","",substr(dir,1,nchar(dir)-1)),"/evaluation_0.2.csv"))
#fwrite(evaluation_high, paste0("../results/sepehri_DMRs/",gsub("^.*\\/","",substr(dir,1,nchar(dir)-1)),"/evaluation_high_0.2.csv"))
#fwrite(evaluation_low, paste0("../results/sepehri_DMRs/",gsub("^.*\\/","",substr(dir,1,nchar(dir)-1)),"/evaluation_low_0.2.csv"))

#fwrite(evaluation, paste0(dir, "evaluation_0.5_normalizations.csv"))

fwrite(evaluation, paste0(dir, "evaluation_0.2.csv"))
fwrite(evaluation_high, paste0(dir, "evaluation_high_0.2.csv"))
fwrite(evaluation_low, paste0(dir, "evaluation_low_0.2.csv"))

## Cluster adjacent probes based on pearson correlation values
## you will need beta values and Illumina 450k manifest file
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)

# return probes of chr
chr.probes <- function(chr){
  ann450k <- fread('data/manifest/450k_manifest.csv')
  ps <- ann450k[ann450k$CHR == chr, 'IlmnID']
  return(ps)
}

# order probes in the genome
order_in_genome <- function(df, probe_names){
  manifest <- fread('data/manifest/450k_manifest.csv', data.table = F)
  manifest <- manifest[!(is.na(manifest$MAPINFO)),]
  manifest <- manifest[,c("IlmnID", "CHR", "MAPINFO")]
  df$IlmnID <- probe_names
  df <- merge(df, manifest, by = c('IlmnID'))
  chrs <- unique(df$CHR)
  sorted_df <- df[FALSE,]
  for (c in chrs){
    temp_df <- df[df$CHR == c,]
    temp_df <- temp_df[order(temp_df$MAPINFO, decreasing = FALSE),]
    sorted_df <- rbind(sorted_df, temp_df)
  }
  rownames(sorted_df) <- sorted_df$IlmnID
  sort_map <- subset(sorted_df, select = c(IlmnID, CHR, MAPINFO))
  sorted_df <- subset(sorted_df, select = -c(IlmnID, CHR, MAPINFO))
  return(list(sort_map, sorted_df))
}

loc_to_probe <- function(start, end, map){
  s_i <- which(map$MAPINFO == start,arr.ind = T)
  e_i <- which(map$MAPINFO == end,arr.ind = T)
  return(list(map$IlmnID[s_i:e_i]))
}

#read beta values of all 1202 samples
data <- fread("data/monocyte_raw_data/GSE56046_b.csv", data.table = F)
names(data) <- data[1,]
data <- data[-1,]
rownames(data) <- data[,1]
data <- data[,-1]

chr1.probes <- chr.probes('1')
data <- data[intersect(chr1.probes$IlmnID, rownames(data)),]
sum(is.na(data))
gc()

# order in genome
sorting <- order_in_genome(data, rownames(data))
sort_map <- sorting[[1]]
sort_df <- sorting[[2]]
sum(rownames(sort_df) != sort_map$IlmnID)

sort_map$correlation <- NA
sort_map$random_correlation <- NA
sort_map$dist <- NA
sort_all <- cbind(sort_map, sort_df)
sort_all_splitted <- split(sort_all, f=sort_all$CHR)
chr_num <- "1"
sort_all_splitted <- sort_all_splitted[chr_num]

# calculate correlation values for each probe with n=window.size neighbors
window.size <- 1

for (m in 1:length(sort_all_splitted)){
  m_df <- sort_all_splitted[[m]][,1:6]
  s_df <- sort_all_splitted[[m]][,7:dim(sort_all_splitted[[m]])[2]]
  for (i in 1:dim(m_df)[1]){
    if (i <= window.size){
      s = -1
      for (j in 1:(i+window.size)){
        s <- s + cor(as.numeric(s_df[i,]), as.numeric(s_df[j,]))
      }
      m_df$correlation[i] <- s/(i+window.size-1)
      next
    }
    
    if (i >= (dim(m_df)[1]-window.size+1)){
      s = -1
      for (j in (i-window.size):(dim(m_df)[1]) ){
        s <- s + cor(as.numeric(s_df[i,]), as.numeric(s_df[j,]))
      }
      m_df$correlation[i] <- s/(dim(m_df)[1]-i+window.size)
      next
    }
    
    s = -1
    for (j in (i-window.size):(i+window.size)){
      s <- s + cor(as.numeric(s_df[i,]), as.numeric(s_df[j,]))
    }
    m_df$correlation[i] <- s/(2*window.size)
  }
  sort_all_splitted[[m]]$correlation <- m_df$correlation
}

# calculate correlatoin values for each probe with n=window.size random probes
for (m in 1:length(sort_all_splitted)){
  m_df <- sort_all_splitted[[m]][,1:6]
  s_df <- sort_all_splitted[[m]][,7:dim(sort_all_splitted[[m]])[2]]
  for (i in 1:dim(m_df)[1]){
    rand_locs <- c(1:dim(m_df)[1])
    rand_locs <- rand_locs[rand_locs != i]
    rand_locs <- sample(rand_locs, 2*window.size)
    s = 0
    for (j in rand_locs){
      s <- s + cor(as.numeric(s_df[i,]), as.numeric(s_df[j,]))
    }
    m_df$random_correlation[i] <- s/(2*window.size)
  }
  sort_all_splitted[[m]]$random_correlation <- m_df$random_correlation
}

for (m in 1:length(sort_all_splitted)){
  m_df <- sort_all_splitted[[m]][,1:6]
  for (i in 2:dim(m_df)[1]){
    m_df$dist[i] <- m_df$MAPINFO[i] - m_df$MAPINFO[i-1]
  }
  sort_all_splitted[[m]]$dist <- m_df$dist
}

sort_all <- bind_rows(sort_all_splitted)
names(sort_all)[4] <- 'correlation'
sum(is.na(sort_all$correlation))

sort_map <- sort_all[,1:6]
sort_df <- sort_all[,7:dim(sort_all)[2]]

fwrite(sort_map, 'results/chr1_sort_map.csv')
fwrite(sort_df, 'results/chr1_sort_df.csv', row.names = T)

sort_map <- fread('results/chr1_sort_map.csv', data.table = F)
sort_df <- fread('results/chr1_sort_df.csv', data.table = F)
rownames(sort_df) <- sort_df[,1]
sort_df <- sort_df[,-1]
names(sort_df) <- sort_df[1,]
sort_df <- sort_df[-1,]

min(sort_map$correlation)
max(sort_map$correlation)
hist(sort_map$correlation)
min(sort_map$random_correlation)
max(sort_map$random_correlation)
hist(sort_map$random_correlation)
quantile(sort_map$correlation)

# compare correlation in neighbors with random probes
plot_df <- sort_map[,c("correlation", "random_correlation")]
plot_df <- melt(plot_df)
names(plot_df) <- c('variable', 'correlation')
mu <- ddply(plot_df, "variable", summarise, grp.mean=mean(correlation))
ggplot(plot_df, aes(x=correlation, fill=variable, color=variable)) +
  geom_histogram(position="identity", alpha=0.5) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=variable),
             linetype="dashed") + theme_bw()

cor_density <- density(plot_df$correlation[plot_df$variable == 'correlation'], from = 0, to = 1)
rand_density <- density(plot_df$correlation[plot_df$variable == 'random_correlation'], from = 0, to = 1)
idx <- (cor_density$y > rand_density$y) &
  (cor_density$x > 0) &
  (cor_density$x < 0.1)
poi <- min(cor_density$x[idx])
poi
max(sort_map$random_correlation)

ggplot(plot_df) +
  geom_density(aes(x = correlation, color = variable)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  scale_color_manual(values = c("tomato", "dodgerblue")) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "values",
       title = "Two density curves") +
  geom_vline(xintercept = poi, linetype = 2, size = 0.3, color = "black") +
  annotate(geom = "text", label = round(poi, 3),
           x = poi - 0.05, y = 0.25, size = 5, angle = 90)

cor_thresh <- 0.1
sort_map$bounds <- as.numeric(sort_map$correlation >= cor_thresh)

plot_df <- sort_map[10:100,]
ggplot(data=plot_df, aes(x=MAPINFO, y=correlation)) +
  geom_point() + geom_line(aes(MAPINFO,bounds))+ theme_bw() +geom_smooth(method = "gam")

# setting filter thresholds based on HMM islands
manifest <- fread('data/manifest/450k_manifest.csv', data.table = F)
manifest <- manifest[!(is.na(manifest$MAPINFO)),]

HMM_islands <- sub('.*:','',manifest$HMM_Island)
HMM_islands_1 <- as.numeric(sub('-.*','',HMM_islands))
HMM_islands_2 <- as.numeric(sub('.*-','',HMM_islands))
HMM_islands_length <- HMM_islands_2 - HMM_islands_1
which(HMM_islands_length < 0)
HMM_islands_length <- HMM_islands_length[!(is.na(HMM_islands_length))]
length(HMM_islands_length)
hist(HMM_islands_length)

plot_df <- data.frame((HMM_islands_length))
names(plot_df) <- c("Length")
ggplot(plot_df, aes(x=Length)) + 
  geom_histogram(binwidth = 50) + 
  geom_vline(aes(xintercept=mean(Length)), color="blue", linetype="dashed", size=1) + theme_bw()

plot_df <- data.frame(HMM_islands_length[HMM_islands_length<10000])
names(plot_df) <- c("Length")
ggplot(plot_df, aes(x=Length)) + 
  geom_histogram(aes(y=..density..), binwidth = 300, colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+ 
  geom_vline(aes(xintercept=mean(HMM_islands_length)), color="blue", linetype="dashed", size=1) + theme_bw()

median(sort_map$dist[!is.na(sort_map$dist)])
mean(HMM_islands_length)
median(HMM_islands_length)
min(HMM_islands_length)
dist_thresh <- median(HMM_islands_length)
length_thresh = min(HMM_islands_length)
dist_thresh <- 702
length_thresh <- 12
remove(HMM_islands)
remove(HMM_islands_1)
remove(HMM_islands_2)
remove(manifest)
gc()

# assign clusters and find location and probes of each cluster
sort_map$clusters <- NA
c = 1
for (i in 1:dim(sort_map)[1]){
  if (sort_map$bounds[i] == 1 & i > 1){
    if (sort_map$bounds[i-1] == 1 & sort_map$dist[i] <= dist_thresh){
      sort_map$clusters[i] <- sort_map$clusters[i-1]
    }
    else{
      sort_map$clusters[i] <- c
      c <- c+1
    }
  }
}

clusters <- sort_map[!is.na(sort_map$clusters),]
write.csv(clusters, 'results/0.1_clusters.csv')
clusters <- read.csv('results/0.1_clusters.csv')
rownames(clusters) <- clusters$X
clusters$X <- NULL

clusters_gpd <- clusters %>% group_by(clusters)
clusters_loc <- data.frame(clusters_gpd %>% dplyr::summarise(start = min(MAPINFO), end = max(MAPINFO), length = end-start))
clusters_loc <- clusters_loc[order(clusters_loc$start),]

selected_clusters_loc <- clusters_loc[clusters_loc$length >= length_thresh,]

manifest <- fread('data/manifest/450k_manifest.csv', data.table = F)
manifest <- manifest[!(is.na(manifest$MAPINFO)),]
rownames(manifest) <- manifest$IlmnID
map <- manifest[manifest$CHR == "1",]
map <- map[order(map$MAPINFO),]
cluster_probes <- c()
cluster_numbers <- c()
for (i in 1:dim(selected_clusters_loc)[1]){
  ps <- loc_to_probe(selected_clusters_loc$start[i],
                     selected_clusters_loc$end[i],
                     map)
  selected_clusters_loc$probes[i] <- ps
  cluster_probes <- c(cluster_probes, ps[[1]])
  selected_clusters_loc$probe_num[i] <- length(selected_clusters_loc$probes[i][[1]])
  cluster_numbers <- c(cluster_numbers, rep(i, selected_clusters_loc$probe_num[i]))
}

selected_clusters_probes <- data.frame(cbind(cluster_numbers, cluster_probes))
dim(selected_clusters_probes)[1]
write.csv(selected_clusters_probes, 'results/chr1_0.1_cluster_probes.csv', row.names = T)

# visualize a DMR
dmr_1 <- data.frame(t(sort_df[selected_clusters_loc[selected_clusters_loc$clusters == 1329, 'probes'][[1]],
                                  sample.int(dim(sort_df)[2], 10)]))
dmr_1$Sample <- rownames(dmr_1)
dmr_1_melted <- melt(dmr_1, id.vars = 'Sample')

ggplot(dmr_1_melted, aes(x = Sample, y = value)) +
  geom_line(aes(color = variable, group = variable)) + geom_point() + theme_bw()

ggplot(dmr_1_melted, aes(x = variable, y = value)) +
  geom_line(aes(color = Sample, group = Sample)) + geom_point() + theme_bw()

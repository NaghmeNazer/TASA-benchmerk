## The code to generate DMRs with seqlm method
## You need genome information downloaded from https://github.com/raivokolde/seqlm
## You need matrix of methylation values (beta here) and phenotypes of interest
## results are saved for further analysis
library(seqlm)
library(data.table)

#Seqlm
load("../data/monocytes/manifest/genome_information.RData")

sample_info <- fread("../data/small_simulated/large_data_info.csv")

dir = "../data/simulated/"
files <- list.files(dir)
files <- files[files != "small_data_info.csv"]
files <- files[grepl("betaQN", files)]
f <- files[1]

for (f in files){
  data <- fread(paste0(dir, f))
  rownames(data) <- data$V1
  data$V1 <- NULL
  print(f)
  print(sum(sample_info$Sample_Name != names(data)))
  segments = seqlm(values = data, genome_information = genome_information, annotation =  sample_info$Tissue)
  dmrs <- data.frame(segments)
  fwrite(dmrs, paste0("../results/seqlm_DMRs/small_CD8/",f))
  rm(list=ls())
  .rs.restartR()
}


### script to sum up evaluation metrics for different tools
library(data.table)

### this function prints best pipeline for each criteria, for example F1, Recall, ...
print_scenario <- function(e_list, scenario = "Best_F1"){
  print(scenario)
  for (i in 1:length(e_list)){
    temp <- e_list[[i]]
    if (scenario == "Best_F1"){
      temp <- temp[order(-temp$F1),]
    }
    if (scenario == "Best_Precision"){
      temp <- temp[order(-temp$Precision),]
    }
    if (scenario == "Best_Recall"){
      temp <- temp[order(-temp$Recall),]
    }
    if (scenario == "Worst_F1"){
      temp <- temp[order(temp$F1),]
    }
    if (scenario == "Worst_Precision"){
      temp <- temp[order(temp$Precision),]
    }
    if (scenario == "Worst_Recall"){
      temp <- temp[order(-temp$Recall),]
    }
    print(paste0("Comparison_Scenario : ", temp[1, "Comparison_Scenario"]))
    print("Best Practice is as follows:")
    print(paste0("Normalization Method = ", temp[1, "Normalization"]))
    print(paste0("Batch Effect Correction = ", temp[1, "BatchEffectCorrected"]))
    print(paste0("DMR Method = ", temp[1, "DMR_Method"]))
    print("*****************************")
  }
}

#evaluation_0.5_dirs <- list.files("../results", "evaluation_0.5.csv", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
evaluation_0.2_dirs <- list.files("../results", "evaluation_0.2.csv", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
f <- evaluation_0.2_dirs[1]
total_eval <- fread(f, data.table = F)

for (i in 2:length(evaluation_0.2_dirs)){
  f <- evaluation_0.2_dirs[i]
  total_eval <- rbind(total_eval, fread(f, data.table = F))
}

table(total_eval$Method)

total_eval$F1 <- (2*total_eval$Precision*total_eval$Recall)/(total_eval$Precision+total_eval$Recall)

total_eval$Tissue <- "CD8"
total_eval$Tissue[grepl("Breast", total_eval$Scenario)] <- "Breast"
total_eval$Tissue[grepl("ns", total_eval$File) | grepl("no_simulation", total_eval$File)] <- "No_Sim"

total_eval$Size <- "Large"
total_eval$Size[grepl("small", total_eval$Scenario)] <- "Small"

total_eval$Threshold <- substr(total_eval$File, 1, 3)
total_eval$Threshold[grepl("ns", total_eval$File) | grepl("no_simulation", total_eval$File)] <- "0"

total_eval$Normalization <- "Raw"
total_eval$Normalization[grepl("BAQN", total_eval$File)] <-"BAQN"

total_eval$BatchEffectCorrected <- "YES"
total_eval$BatchEffectCorrected[grepl("not",total_eval$File) | grepl("bnc", total_eval$File)] <- "NO"

total_eval$Comparison <- paste0(total_eval$Size, "_", total_eval$Tissue, "_", total_eval$Threshold)
table(total_eval$Comparison)

total_eval <- total_eval[,c("Comparison", "Normalization", "BatchEffectCorrected", "Method",
                            "TP", "FP", "FN", "Precision", "Recall", "F1")]

unique_scenario <- unique(total_eval$Comparison)
comparison_list <- list()
for (s in unique_scenario){
  temp <- total_eval[total_eval$Comparison == s,]
  if (grepl("No_Sim", s)){
    temp <- temp[order(temp$FP),]
  }
  else{
    temp <- temp[order(-temp$F1),]
  }
  temp[is.na(temp)] <- 0
  names(temp) <- c("Comparison_Scenario", "Normalization", "BatchEffectCorrected", "DMR_Method",
                    "TP", "FP", "FN", "Precision", "Recall", "F1")
  fwrite(temp,paste0("../results/evaluation_0.2/", s, ".csv"))
  comparison_list[[s]] <- temp
}

print_scenario(comparison_list, "Best_F1")
print_scenario(comparison_list, "Best_Precision")
print_scenario(comparison_list, "Best_Recall")
print_scenario(comparison_list, "Worst_F1")
print_scenario(comparison_list, "Worst_Precision")
print_scenario(comparison_list, "Worst_Recall")



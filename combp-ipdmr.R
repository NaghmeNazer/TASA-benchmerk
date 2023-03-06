## the code to generate DMRs with combp and ipdmr methods
## both need a matrix of p values for each probe, here we used the dmps generated with dmpfinder function of minfi
## the results of saved for further analysis
library(ENmix)
library(data.table)

manifest <- fread("/media/mehrmohammadi_hdd/nazer/data/monocytes/manifest/450k_manifest.csv", data.table = F)
rownames(manifest) <- manifest$IlmnID
idmr_probes <- fread("/media/mehrmohammadi_hdd/nazer/data/monocytes/control_probes/iDMR_probes.csv", data.table = F)
idmr_probes <- as.character(idmr_probes[-1,-1])

dir = "/media/mehrmohammadi_hdd/sepehri/sepehri/sepehri/simulation/minfi_results/CD8_small//"
#dir = "/home/sepehri/simulation/minfi_results/CD8_large/"
files <- list.files(dir)
#files <- files[grepl("ns", files)]
files
res_dir = "small_CD8"

for (f in files){
  dmps <- fread(paste0(dir, f), data.table = F)
  names(dmps) <- c("probe", "intercept", "f", "p", "q")
  dmps <- dmps[!(dmps$probe %in% idmr_probes),]
  dmps$chr <- manifest[dmps$probe, "CHR"]
  dmps$start <- manifest[dmps$probe, "MAPINFO"]
  dmps$end <- dmps$start + 1
  dmps <- dmps[,c("chr", "start", "end", "p", "probe")]
  
  ## run combp
  combp(data=dmps,seed=0.05, region_plot = F, mht_plot = F)
  combp_dmrs <- fread("resu_combp.csv", data.table = F)
  fwrite(combp_dmrs, paste0("../results/combp_DMRs/",res_dir,"/",f))
  
  ## run ipdmr
  ipdmr(data=dmps,seed=0.05, region_plot = F, mht_plot = F)
  ipdmr_dmrs <- fread("resu_ipdmr.csv", data.table = F)
  fwrite(ipdmr_dmrs, paste0("../results/ipdmr_DMRs/",res_dir,"/",f))
}


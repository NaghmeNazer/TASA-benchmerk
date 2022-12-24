library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(wateRmelon)

#####
#find infinium types of probes to separate probes for normalization
get_infinium_type <- function(cgs){
  m <- fread('data/manifest/450k_manifest.csv', data.table = F)
  rownames(m) <- m$IlmnID
  m <- m[cgs, c("IlmnID", "Infinium_Design_Type")]
  m$Infinium_Design_Type[m$Infinium_Design_Type == 'I'] <- 1
  m$Infinium_Design_Type[m$Infinium_Design_Type == 'II'] <- 2
  return(m)
}

#BMIQ normalization
modified_BMIQ <- function(beta.v,design.v,nL=3,doH=TRUE,nfit=10000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1,pri=TRUE){
  ##LS>
  if(!library(RPMM, logical.return=TRUE, quietly=TRUE)){
    stop( 'need RPMM package')
  }
  good     <- !is.na(beta.v)
  out      <- beta.v
  beta.v   <- beta.v[good]
  design.v <- design.v[good]
  print    <- function(x){ if(pri)base::print(x)}
  ##LS<
  type1.idx <- which(design.v==1);
  type2.idx <- which(design.v==2);
  
  beta1.v <- beta.v[type1.idx];
  beta2.v <- beta.v[type2.idx];
  
  ### check if there are exact 0's or 1's. If so, regularise using minimum positive and maximum below 1 values.
  if(min(beta1.v)==0){
    beta1.v[beta1.v==0] <- min(setdiff(beta1.v,0));
  }
  if(min(beta2.v)==0){
    beta2.v[beta2.v==0] <- min(setdiff(beta2.v,0));
  }
  if(max(beta1.v)==1){
    beta1.v[beta1.v==1] <- max(setdiff(beta1.v,1));
  }
  if(max(beta2.v)==1){
    beta2.v[beta2.v==1] <- max(setdiff(beta2.v,1));
  }
  
  ### estimate initial weight matrix from type1 distribution
  w0.m <- matrix(0,nrow=length(beta1.v),ncol=nL);
  w0.m[which(beta1.v <= th1.v[1]),1] <- 1;
  w0.m[intersect(which(beta1.v > th1.v[1]),which(beta1.v <= th1.v[2])),2] <- 1;
  w0.m[which(beta1.v > th1.v[2]),3] <- 1;
  
  ### fit type1
  print("Fitting EM beta mixture to type1 probes");
  rand.idx <- sample(1:length(beta1.v),nfit,replace=FALSE)
  #rand.idx <- 1:length(beta1.v)
  em1.o <- blc(matrix(beta1.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol,verbose=pri); ##LS 
  subsetclass1.v <- apply(em1.o$w,1,which.max);
  subsetth1.v <- c(mean(max(beta1.v[rand.idx[subsetclass1.v==1]]),min(beta1.v[rand.idx[subsetclass1.v==2]])),mean(max(beta1.v[rand.idx[subsetclass1.v==2]]),min(beta1.v[rand.idx[subsetclass1.v==3]])));
  subsetth1.v <- th1.v
  class1.v <- rep(2,length(beta1.v));
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1;
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3;
  nth1.v <- subsetth1.v;
  print("Done");
  
  ### generate plot from estimated mixture
  if(plots){
    print("Check");
    tmpL.v <- as.vector(rmultinom(1:nL,length(beta1.v),prob=em1.o$eta));
    tmpB.v <- vector();
    for(l in 1:nL){
      tmpB.v <- c(tmpB.v,rbeta(tmpL.v[l],em1.o$a[l,1],em1.o$b[l,1]));
    }
    
    pdf(paste("Type1fit-",sampleID,".pdf",sep=""),width=6,height=4);
    plot(density(beta1.v));
    d.o <- density(tmpB.v);
    points(d.o$x,d.o$y,col="green",type="l")
    legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
    dev.off();
  }
  
  
  
  ### Estimate Modes 
  d1U.o <- density(beta1.v[class1.v==1])
  d1M.o <- density(beta1.v[class1.v==3])
  mod1U <- d1U.o$x[which.max(d1U.o$y)]
  mod1M <- d1M.o$x[which.max(d1M.o$y)]
  d2U.o <- density(beta2.v[which(beta2.v<0.4)]);
  d2M.o <- density(beta2.v[which(beta2.v>0.6)]);
  mod2U <- d2U.o$x[which.max(d2U.o$y)]
  mod2M <- d2M.o$x[which.max(d2M.o$y)]
  
  
  ### now deal with type2 fit
  th2.v <- vector();
  th2.v[1] <- nth1.v[1] + (mod2U-mod1U);
  th2.v[2] <- nth1.v[2] + (mod2M-mod1M);
  
  ### estimate initial weight matrix 
  w0.m <- matrix(0,nrow=length(beta2.v),ncol=nL);
  w0.m[which(beta2.v <= th2.v[1]),1] <- 1;
  w0.m[intersect(which(beta2.v > th2.v[1]),which(beta2.v <= th2.v[2])),2] <- 1;
  w0.m[which(beta2.v > th2.v[2]),3] <- 1;
  
  print("Fitting EM beta mixture to type2 probes");
  rand.idx <- sample(1:length(beta1.v),nfit,replace=FALSE)
  em2.o <- blc(matrix(beta2.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol,verbose=pri); ##LS
  print("Done");
  
  ### for type II probes assign to state (unmethylated, hemi or full methylation)
  subsetclass2.v <- apply(em2.o$w,1,which.max);
  #subsetth2.v <- c(mean(max(beta2.v[rand.idx[subsetclass2.v==1]]),min(beta2.v[rand.idx[subsetclass2.v==2]])),mean(max(beta2.v[rand.idx[subsetclass2.v==2]]),min(beta2.v[rand.idx[subsetclass2.v==3]])));
  subsetth2.v <- th2.v
  class2.v <- rep(2,length(beta2.v));
  class2.v[which(beta2.v < subsetth2.v[1])] <- 1;
  class2.v[which(beta2.v > subsetth2.v[2])] <- 3;
  
  
  ### generate plot
  if(plots){
    tmpL.v <- as.vector(rmultinom(1:nL,length(beta2.v),prob=em2.o$eta));
    tmpB.v <- vector();
    for(lt in 1:nL){
      tmpB.v <- c(tmpB.v,rbeta(tmpL.v[lt],em2.o$a[lt,1],em2.o$b[lt,1]));
    }
    pdf(paste("Type2fit-",sampleID,".pdf",sep=""),width=6,height=4);
    plot(density(beta2.v));
    d.o <- plot(density(tmpB.v));
    points(d.o$x,d.o$y,col="green",type="l")
    legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
    dev.off();
  }
  
  classAV1.v <- vector();classAV2.v <- vector();
  for(l in 1:nL){
    classAV1.v[l] <-  em1.o$mu[l,1];
    classAV2.v[l] <-  em2.o$mu[l,1];
  }
  
  ### start normalising type2 probes
  print("Start normalising type 2 probes");
  nbeta2.v <- beta2.v;
  ### select U probes
  lt <- 1;
  selU.idx <- which(class2.v==lt);
  selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])];
  selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])];
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selUR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
  nbeta2.v[selUR.idx] <- q.v;
  p.v <- pbeta(beta2.v[selUL.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=TRUE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=TRUE);
  nbeta2.v[selUL.idx] <- q.v;
  
  ### select M probes
  lt <- 3;
  selM.idx <- which(class2.v==lt);
  selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])];
  selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])];
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selMR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
  nbeta2.v[selMR.idx] <- q.v;
  
  
  if(doH){ ### if TRUE also correct type2 hemimethylated probes
    ### select H probes and include ML probes (left ML tail is not well described by a beta-distribution).
    lt <- 2;
    selH.idx <- c(which(class2.v==lt),selML.idx);
    minH <- min(beta2.v[selH.idx])
    maxH <- max(beta2.v[selH.idx])
    deltaH <- maxH - minH;
    #### need to do some patching
    deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
    deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])
    
    ## new maximum of H probes should be
    nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM;
    ## new minimum of H probes should be
    nminH <- max(nbeta2.v[selU.idx]) + deltaUH;
    ndeltaH <- nmaxH - nminH;
    
    ### perform conformal transformation (shift+dilation)
    ## new_beta_H(i) = a + hf*(beta_H(i)-minH);
    hf <- ndeltaH/deltaH ;
    ### fix lower point first
    nbeta2.v[selH.idx] <- nminH + hf*(beta2.v[selH.idx]-minH);
    
  }
  
  pnbeta.v <- beta.v;
  pnbeta.v[type1.idx] <- beta1.v;
  pnbeta.v[type2.idx] <- nbeta2.v;
  
  ### generate final plot to check normalisation
  if(plots){
    print("Generating final plot");
    d1.o <- density(beta1.v);
    d2.o <- density(beta2.v);
    d2n.o <- density(nbeta2.v);
    ymax <- max(d2.o$y,d1.o$y,d2n.o$y);
    pdf(paste("CheckBMIQ-",sampleID,".pdf",sep=""),width=6,height=4)
    plot(density(beta2.v),type="l",ylim=c(0,ymax),xlim=c(0,1));
    points(d1.o$x,d1.o$y,col="red",type="l");
    points(d2n.o$x,d2n.o$y,col="blue",type="l");
    legend(x=0.5,y=ymax,legend=c("type1","type2","type2-BMIQ"),bty="n",fill=c("red","black","blue"));
    dev.off();
  }
  
  print(paste("Finished for sample ",sampleID,sep=""));
  
  ##LS>
  out[good] <- pnbeta.v
  pnbeta.v  <- out
  ##LS<
  
  return(list(nbeta=pnbeta.v,class1=class1.v,class2=class2.v,av1=classAV1.v,av2=classAV2.v,hf=hf,th1=nth1.v,th2=th2.v));
}

# BAQN (Between-array Qunatile Normalization)
baqn <- function(df, infinium_types){
  selected_df_BAQN <- df
  selected_df_BAQN_1 <-  selected_df_BAQN[infinium_types[infinium_types$Infinium_Design_Type == 1, 'IlmnID'],]
  selected_df_BAQN_2 <-  selected_df_BAQN[infinium_types[infinium_types$Infinium_Design_Type == 2, 'IlmnID'],]
  selected_df_BAQN_1 <- normalizeQuantiles(selected_df_BAQN_1)
  selected_df_BAQN_2 <- normalizeQuantiles(selected_df_BAQN_2)
  selected_df_BAQN <- rbind(selected_df_BAQN_1, selected_df_BAQN_2)
  return(selected_df_BAQN)
}
#####
## reading control probes to evaluate algorithms' efficiency
control.probes <- fread('data/control_probes/control_probes.csv', data.table = F)
control.probes <- control.probes[-1,2]
idmr.probes <- fread('data/control_probes/iDMR_probes.csv', data.table = F)
idmr.probes <- idmr.probes[-1,2]
idmr0 <- fread('data/monocyte_raw_data/GSE56046_b.csv', data.table = F)
names(idmr0) <- idmr0[1,]
idmr0 <- idmr0[-1,]
rownames(idmr0) <- idmr0[,1]
idmr0 <- idmr0[,-1]
idmr0 <- idmr0[idmr.probes,]
idmr0 <- idmr0[complete.cases(idmr0),]
sum(is.na(idmr0))

idmr1 <- fread('data/control_probes/iDMR_probes_GSE120610.csv', data.table = F)
rownames(idmr1) <- idmr1[,1]
idmr1 <- idmr1[,-1]
sum(is.na(idmr1))
names(idmr1) <- sub(".Methylated.Signal", "", names(idmr1))

idmr2 <- fread('data/control_probes/iDMR_probes_GSE131989.csv', data.table = F)
rownames(idmr2) <- idmr2[,1]
idmr2 <- idmr2[,-1]
sum(is.na(idmr2))

idmr3 <- fread('data/control_probes/iDMR_probes_GSE134429.csv', data.table = F)
rownames(idmr3) <- idmr3[,1]
idmr3 <- idmr3[,-1]
sum(is.na(idmr3))

idmr4 <- fread('data/control_probes/iDMR_probes_GSE184269.csv', data.table = F)
rownames(idmr4) <- idmr4[,1]
idmr4 <- idmr4[,-1]
sum(is.na(idmr4))

idmr.intersect.probes <- Reduce(intersect, list(rownames(idmr0),
                                                rownames(idmr1),
                                                rownames(idmr2),
                                                rownames(idmr3),
                                                rownames(idmr4)))
#####
## read data -- choose appropriate data based on simulation scenario
data_name <- '../data/simulated/0.4_large_data_raw_batch_corrected.csv'
data <- fread(data_name, data.table = F)
rownames(data) <- data[,1]
data <- data[,-1]
names(data) <- gsub(".Methylated.Signal","", names(data))

## filter idmr data
idmr_total <- cbind(idmr0[idmr.intersect.probes,],
                    idmr1[idmr.intersect.probes,],
                    idmr2[idmr.intersect.probes,],
                    idmr3[idmr.intersect.probes,],
                    idmr4[idmr.intersect.probes,])
sum(is.na(idmr_total))
idmr_total <- idmr_total[,names(data)]
setdiff(names(data), names(idmr_total))

## Normalization
selected_df <- rbind(data, idmr_total)
infinium_types <- get_infinium_type(rownames(selected_df))

# Perform BMIQ normalization
selected_df_BMIQ <- selected_df
selected_df_BMIQ[,] <- NA
for (s in names(selected_df)){
  selected_df_BMIQ[,s] <- modified_BMIQ(selected_df[,s], infinium_types$Infinium_Design_Type, nfit = 10000, plots = F)[["nbeta"]]
}

# Perform betaQN normalization
selected_df_betaqn <- betaqn(selected_df)

# Perform BAQN normalization
selected_df_BAQN <- baqn(selected_df, infinium_types)

## Compute evaluation metrics
dmrse_row(as.matrix(selected_df))
dmrse_row(as.matrix(selected_df_BMIQ))
dmrse_row(as.matrix(selected_df_betaqn))
dmrse_row(as.matrix(selected_df_BAQN))

normalization_sds <- data.frame(matrix(,nrow = dim(selected_df)[1], ncol = 4))
names(normalization_sds) <- c('Raw', 'BMIQ', 'BetaQN', 'BAQN')
rownames(normalization_sds) <- rownames(selected_df)
normalization_sds$Raw <- rowSds(as.matrix(selected_df))
normalization_sds$BMIQ <- rowSds(as.matrix(selected_df_BMIQ))
normalization_sds$BetaQN <- rowSds(as.matrix(selected_df_betaqn))
normalization_sds$BAQN <- rowSds(as.matrix(selected_df_BAQN))
normalization_sds$probe_type <- infinium_types$Infinium_Design_Type
plot_df <- normalization_sds
plot_df <- melt(plot_df)
plot_df$category <- paste0(plot_df$variable, '-Type', plot_df$probe_type)

SDMedian <- summarise(group_by(plot_df, category), MD = median(value))
SDMedian <- summarise(group_by(plot_df, variable), MD = median(value))
SDMedian

## write data
fwrite(selected_df_BAQN, gsub('raw', 'BAQN', data_name), col.names = T, row.names = T)

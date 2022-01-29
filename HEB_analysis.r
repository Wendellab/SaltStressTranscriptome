###############################################
###########Homoeolog expresson bias############
###############################################
###load packages
library(DESeq2)
library(tidyverse)
##read in data
load("~/Documents/counts.sub.RData")
# make coldata
species <- c(rep("A2",6),rep("D5",6),rep("AD1",12),rep("AD4",12))
condition <- c(rep(c(rep("Control",3),rep("Salt",3)),2),rep(c(rep("Control",6),rep("Salt",6)),2))
rep <- c(rep(1:3,4),rep(c(1,1,2,2,3,3),4))
colsub <- data.frame(row.names=names(counts.sub),species,condition,rep)
colsub$homoeolog <- c(rep("A",6),rep("D",6),rep(c("A","D"),12))
colsub$ploidy <- c(rep("Diploid",12),rep("Tetraploid",24))
colsub$sample <- paste(colsub$species,colsub$condition,sep = ".")
colsub$base <- paste(colsub$sample,colsub$homoeolog,sep=".")
all(colnames(counts.sub)==rownames(colsub))

# remove D5.ck.2, which was shown to be an outliner
counts.sub29 <- counts.sub[,-8];colsub29 <- colsub[-8,]

# make data for AD1 and AD4 separately
counts.sub29AD1 <- counts.sub29[,-c(24:35)]
counts.sub29AD4 <- counts.sub29[,-c(12:23)]
colsub29AD1 <- colsub29[-c(24:35),]
colsub29AD4 <- colsub29[-c(12:23),]
colsub29AD1$ploidy <- factor(colsub29AD1$ploidy,levels = c("Tetraploid","Diploid"))
colsub29AD4$ploidy <- factor(colsub29AD4$ploidy,levels = c("Tetraploid","Diploid"))

### AD1 analysis
dds.AD1.hpc <- DESeqDataSetFromMatrix(countData = counts.sub29AD1,colData = colsub29AD1,design = ~homoeolog + ploidy + condition + homoeolog:ploidy + homoeolog:condition + ploidy:condition + homoeolog:ploidy:condition)
dds.AD1.hpc$ploidy #check the base levels, make sure polyploid is the reference level
des.AD1.hpc <- DESeq(dds.AD1.hpc,test = "Wald")
contrasts_ad1 <- sapply(resultsNames(des.AD1.hpc),list)

# homoeolog*ploidy*condition, novel HEB in response to stress in polyploid
res_hpc.AD1 <- data.frame(results(des.AD1.hpc,contrast = contrasts_ad1[8],pAdjustMethod = "BH"))
res_hpc.AD1[,gsub("padj","sig",colnames(res_hpc.AD1)[grep("padj", colnames(res_hpc.AD1))])]<-ifelse(res_hpc.AD1$padj>0.05 | is.na(res_hpc.AD1$padj),"","sig")
# homoeolog*condition, HEB in response to stress inherited diploids
res_hc.AD1 <- data.frame(results(des.AD1.hpc,contrast = contrasts_ad1[6],pAdjustMethod = "BH"))
res_hc.AD1[,gsub("padj","sig",colnames(res_hc.AD1)[grep("padj", colnames(res_hc.AD1))])]<-ifelse(res_hc.AD1$padj>0.05 | is.na(res_hc.AD1$padj),"","sig")
# homoeolog*ploidy, novel HEB not respond to stress in polyploid
res_hp.AD1 <- data.frame(results(des.AD1.hpc,contrast = contrasts_ad1[5],pAdjustMethod = "BH"))
res_hp.AD1[,gsub("padj","sig",colnames(res_hp.AD1)[grep("padj", colnames(res_hp.AD1))])]<-ifelse(res_hp.AD1$padj>0.05 | is.na(res_hp.AD1$padj),"","sig")
# homoeolog, HEB not respond to stress and inherited from diploids
res_h.AD1 <- data.frame(results(des.AD1.hpc,contrast = contrasts_ad1[2],pAdjustMethod = "BH"))
res_h.AD1[,gsub("padj","sig",colnames(res_h.AD1)[grep("padj", colnames(res_h.AD1))])]<-ifelse(res_h.AD1$padj>0.05 | is.na(res_h.AD1$padj),"","sig")

colnames(res_hpc.AD1) <- paste("homoe.ploid.trt",colnames(res_hpc.AD1),sep = "_")
colnames(res_hc.AD1) <- paste("homoe.trt",colnames(res_hc.AD1),sep = "_")
colnames(res_hp.AD1) <- paste("homoe.ploid",colnames(res_hp.AD1),sep = "_")
colnames(res_h.AD1) <- paste("homoe",colnames(res_h.AD1),sep = "_")

AD1_all <- cbind(res_hpc.AD1,res_hc.AD1,res_hp.AD1,res_h.AD1)
HEBtest.AD1 <- AD1_all[complete.cases(AD1_all),]

## assign significant categories
HEBtest.AD1$both_trt_sig <- with(HEBtest.AD1,ifelse(homoe.ploid.trt_padj <= 0.05 & homoe.trt_padj <= 0.05,"both_trt",0))
HEBtest.AD1$h_p_trt_sig <- with(HEBtest.AD1,ifelse(homoe.ploid.trt_padj <= 0.05 & homoe.trt_padj > 0.05, "novel_trt",0))
HEBtest.AD1$h_trt_sig <- with(HEBtest.AD1,ifelse(both_trt_sig==0 & h_p_trt_sig==0 & homoe.trt_padj <= 0.05, "inher_trt",0))
HEBtest.AD1$both_sig <- with(HEBtest.AD1,ifelse(both_trt_sig==0 & h_p_trt_sig==0 & h_trt_sig==0 & homoe.ploid_padj < 0.05 & homoe_padj < 0.05, "both",0))
HEBtest.AD1$h_p_sig <- with(HEBtest.AD1,ifelse(both_trt_sig==0 & h_p_trt_sig==0 & both_sig==0 & homoe.ploid_padj < 0.05, "novel",0))
HEBtest.AD1$h_sig <- with(HEBtest.AD1,ifelse(both_trt_sig==0 & h_trt_sig==0 & both_sig==0 & homoe_padj < 0.05, "inher",0))

HEBtest.AD1$sig <- with(HEBtest.AD1,paste(both_trt_sig,h_p_trt_sig,h_trt_sig,both_sig,h_p_sig,h_sig, sep=""))
HEBtest.AD1$sig <- gsub("0","",HEBtest.AD1$sig)
HEBtest.AD1$sig[HEBtest.AD1$sig==""] <- "nonsig"
table(HEBtest.AD1$sig)
# both       both_trt          inher      inher_trt inher_trtnovel         nonsig 
#6703            290           4738            166            223           5956 
#novel      novel_trt novel_trtinher 
#4787           1661           1783 

### AD4 analysis
dds.AD4.hpc <- DESeqDataSetFromMatrix(countData = counts.sub29AD4,colData = colsub29AD4,design = ~ homoeolog + ploidy + condition + homoeolog:ploidy + homoeolog:condition + ploidy:condition + homoeolog:ploidy:condition)
dds.AD4.hpc$ploidy
des.AD4.hpc <- DESeq(dds.AD4.hpc,test="Wald")
contrasts_ad4 <- sapply(resultsNames(des.AD4.hpc),list)

# homoeolog*ploidy*condition
res_hpc.AD4 <- data.frame(results(des.AD4.hpc,contrast = contrasts_ad4[8],pAdjustMethod = "BH"))
res_hpc.AD4[,gsub("padj","sig",colnames(res_hpc.AD4)[grep("padj", colnames(res_hpc.AD4))])]<-ifelse(res_hpc.AD4$padj>0.05 | is.na(res_hpc.AD4$padj),"","sig")
# homoeolog*condition
res_hc.AD4 <- data.frame(results(des.AD4.hpc,contrast = contrasts_ad4[6],pAdjustMethod = "BH"))
res_hc.AD4[,gsub("padj","sig",colnames(res_hc.AD4)[grep("padj", colnames(res_hc.AD4))])] <- ifelse(res_hc.AD4$padj>0.05 | is.na(res_hc.AD4$padj),"","sig")
# homoeolog*ploidy
res_hp.AD4 <- data.frame(results(des.AD4.hpc,contrast = contrasts_ad4[5],pAdjustMethod = "BH")) 
res_hp.AD4[,gsub("padj","sig",colnames(res_hp.AD4)[grep("padj",colnames(res_hp.AD4))])] <- ifelse(res_hp.AD4$padj>0.05 | is.na(res_hp.AD4$padj),"","sig")
# homoeolog
res_h.AD4 <- data.frame(results(des.AD4.hpc,contrast = contrasts_ad4[2],pAdjustMethod = "BH"))
res_h.AD4[,gsub("padj","sig",colnames(res_h.AD4)[grep("padj",colnames(res_h.AD4))])] <- ifelse(res_h.AD4$padj>0.05 | is.na(res_h.AD4$padj),"","sig")

colnames(res_hpc.AD4) <- paste("homoe.ploid.trt",colnames(res_hpc.AD4),sep = "_")
colnames(res_hc.AD4) <- paste("homoe.trt",colnames(res_hc.AD4),sep = "_")
colnames(res_hp.AD4) <- paste("homoe.ploid",colnames(res_hp.AD4),sep = "_")
colnames(res_h.AD4) <- paste("homoe",colnames(res_h.AD4),sep = "_")

AD4_all <- cbind(res_hpc.AD4,res_hc.AD4,res_hp.AD4,res_h.AD4)
HEBtest.AD4 <- AD4_all[complete.cases(AD4_all),]

## assign significant categories
HEBtest.AD4$both_trt_sig <- with(HEBtest.AD4,ifelse(homoe.ploid.trt_padj <= 0.05 & homoe.trt_padj <= 0.05,"both_trt",0))
HEBtest.AD4$h_p_trt_sig <- with(HEBtest.AD4,ifelse(homoe.ploid.trt_padj <= 0.05 & homoe.trt_padj > 0.05, "novel_trt",0))
HEBtest.AD4$h_trt_sig <- with(HEBtest.AD4,ifelse(both_trt_sig==0 & h_p_trt_sig==0 & homoe.trt_padj <= 0.05, "inher_trt",0))
HEBtest.AD4$both_sig <- with(HEBtest.AD4,ifelse(both_trt_sig==0 & h_p_trt_sig==0 & h_trt_sig==0 & homoe.ploid_padj < 0.05 & homoe_padj < 0.05, "both",0))
HEBtest.AD4$h_p_sig <- with(HEBtest.AD4,ifelse(both_trt_sig==0 & h_p_trt_sig==0 & both_sig==0 & homoe.ploid_padj < 0.05, "novel",0))
HEBtest.AD4$h_sig <- with(HEBtest.AD4,ifelse(both_trt_sig==0 & h_trt_sig==0 & both_sig==0 & homoe_padj < 0.05, "inher",0))

HEBtest.AD4$sig <- with(HEBtest.AD4,paste(both_trt_sig,h_p_trt_sig,h_trt_sig,both_sig,h_p_sig,h_sig, sep=""))
HEBtest.AD4$sig <- gsub("0","",HEBtest.AD4$sig)
HEBtest.AD4$sig[HEBtest.AD4$sig==""] <- "nonsig"
table(HEBtest.AD4$sig)
#both       both_trt          inher      inher_trt inher_trtnovel         nonsig 
#6381            261           4551            137            184           6114 
#novel      novel_trt novel_trtinher 
#4843           1802           1996 

save(HEBtest.AD1,HEBtest.AD4,file="HEB_test.RData")
###

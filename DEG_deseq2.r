##import counts data
counttable <- read.table("salt.counts",header=TRUE, sep="\t"); dim(counttable)
head(counttable);
names(counttable)
names(counttable) <- gsub("sort.bam","sort.T.bam",names(counttable))
names(counttable) <- gsub(".sort|.bam|_5.22|_5.23|_6.11","",names(counttable))
dim(saltcount <- counttable[,-c(2:6)])
rownames(saltcount) <- saltcount[,1]

#get subset data for each species
A2count <- saltcount[,grep("A2",names(saltcount))];dim(A2count)
D5count <- saltcount[,grep("D5",names(saltcount))];dim(D5count)
TM1count <- saltcount[,grep("TM1",names(saltcount))];dim(TM1count)
PS7count <- saltcount[,grep("PS7",names(saltcount))];dim(PS7count)
AD4count <- saltcount[,grep("AD4",names(saltcount))];dim(AD4count)
names(A2count)
#combine data from "web" in A2 and PS7
A2ck1 <- A2count[,1:5]+A2count[,6:10]
A2ck3 <- A2count[,16:20]+A2count[,21:25]
A2salt2 <- A2count[,31:35]+A2count[,36:40]
A2salt3 <- A2count[,41:45]+A2count[,46:50]
A2count_c <- cbind(A2ck1,A2count[,grep("CK2",names(A2count))],A2ck3,A2count[,grep("Salt1",names(A2count))],A2salt2,A2salt3);dim(A2count_c)

PS7salt3 <- PS7count[,26:30]+PS7count[,31:35]
PS7count_c <- cbind(PS7count[,-grep("Salt3",names(PS7count))],PS7salt3)

#get the new saltcount data
counts <- cbind(A2count_c,D5count,TM1count,PS7count_c,AD4count);dim(counts)

#pull out the Total counts values
counts.T <- counts[,grep(".T$",names(counts))];dim(counts.T)
head(counts.T)
colSums(counts.T)
save(counts.T, file = "counts.T.RData")

#pull out total counts of A2 and D5
diploid.T <- counts[,grep("A2.*T$|D5.*T$",names(counts))];dim(diploid.T)
#A and D counts of polyploids
TM1.AD <- counts[,grep("TM1.*A$|TM1.*D$",names(counts))];dim(TM1.AD)
PS7.AD <- counts[,grep("PS7.*A$|PS7.*D$",names(counts))];dim(PS7.AD)
AD4.AD <- counts[,grep("AD4.*A$|AD4.*D$",names(counts))];dim(AD4.AD)

#get the count table for analysis, ADs subgenome
counts.sub <- cbind(diploid.T, TM1.AD, PS7.AD, AD4.AD)
names(counts.sub) <- gsub("[.]T","",names(counts.sub))
save(counts.sub, file = "counts.sub.RData")

#At+Dt for allopolyploids
TM1.Total <- TM1.AD[,c(seq(1,11,by=2))]+TM1.AD[,c(seq(2,12,by=2))]
PS7.Total <- PS7.AD[,c(seq(1,11,by=2))]+PS7.AD[,c(seq(2,12,by=2))]
AD4.Total <- AD4.AD[,c(seq(1,11,by=2))]+AD4.AD[,c(seq(2,12,by=2))]
counts.total <- cbind(diploid.T, TM1.Total, PS7.Total, AD4.Total)
names(counts.total) <- gsub("[.]T|[.]A","",names(counts.total))
save(counts.total, file = "counts.total.RData")

##load packages
library(DESeq2)
library(ggplot2)
library(scales)
library(ape)
library(pheatmap)
library( RColorBrewer)

#make a coldata
species <- c(rep("A2",6),rep("D5",6),rep("TM1",6),rep("PS7",6),rep("AD4",6))
condition <- rep(c(rep("Control",3),rep("Salt",3)),5)
rep <- rep(1:3,10)
coltotal <- data.frame(row.names=names(counts.total),species,condition,rep)
coltotal$sample <- paste(coltotal$species,coltotal$condition,sep = ".")
all(colnames(counts.total)==rownames(coltotal))
coltotal$condition <- relevel(coltotal$condition,ref = "Control")
write.table(coltotal, file="coldatatotal.txt",row.names = T)

#Build DESeq data
dds <- DESeqDataSetFromMatrix(countData = counts.total, colData = coltotal, design = ~ species + condition + species:condition)
#salt.dds <- salt.dds[rowSums(counts(salt.dds))/30 >= 1,]
#counts.total.test <- counts.total[!(rowSums(counts.total) %in% 0),] #still has 36453 genes
counts.1 <- counts.total[rowSums(counts.total[,grep("A2_CK",names(counts.total))])>=3 |rowSums(counts.total[,grep("A2_Salt",names(counts.total))])>=3 | rowSums(counts.total[,grep("D5_CK",names(counts.total))])>=2 | rowSums(counts.total[,grep("D5_Salt",names(counts.total))])>=3 | rowSums(counts.total[,grep("TM1_CK",names(counts.total))])>=3 | rowSums(counts.total[,grep("TM1_Salt",names(counts.total))])>=3 | rowSums(counts.total[,grep("PS7_CK",names(counts.total))]>=3 | rowSums(counts.total[,grep("PS7_Salt",names(counts.total))])>=3 | rowSums(counts.total[,grep("AD4_CK",names(counts.total))])>=3 | rowSums(counts.total[,grep("PS7_Salt",names(counts.total))])>=3),]
dim(counts.1)

#rlog normalization
rld <- rlog(dds, blind = FALSE)
count.rld <- as.data.frame(assay(rld))

mytheme <- theme_bw()+
  #theme_classic()+
  # scale_color_manual(values = mi, guide = guide_legend(title = NULL))+
  # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
  theme(
    
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    
    plot.title = element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =element_text(size = 18,face = "bold",colour = "black"),
    axis.title.x =element_text(size = 20,face = "bold",colour = "black"),
    axis.text = element_text(size = 16,face = "bold"),
    axis.text.x = element_text(colour = "black",size = 20),
    axis.text.y = element_text(colour = "black",size = 20),
    legend.text = element_text(size = 12,face = "bold"))+
  #theme(legend.position = c(0.1,0.2))+
  
  theme(strip.text.x = element_text(size=15, angle=0),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="blue", fill="#CCCCFF"))

#PCA plot
pca=prcomp(t(count.rld))
dat = as.data.frame(pca$x)
proportion<-summary(pca)$importance[2,1:2]*100
proportion<-paste0(names(proportion)," (", proportion, "%)")
#save picture
tiff(file="PCA_rld.tiff",width = 800,height = 800,units = "px")
ggplot(aes(PC1, PC2, color=species, shape=condition),data=dat) + geom_point(size=5) +xlab(proportion[1]) + ylab(proportion[2])+mytheme
dev.off()

#heatmap
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
head(rld_cor)
tiff(file="correlation_rld.tiff",width = 800,height = 800)
pheatmap(rld_cor)
dev.off()

#remove D5_Ck2 which shown to be an outliner
coltotal29 <- coltotal[-8,]
dds <- DESeqDataSetFromMatrix(countData = counts.total[,-8], colData = coltotal29, design = ~ species + condition + species:condition)
dds$condition
dds$species

dds <- estimateSizeFactors(dds)
ddstotal <- DESeq(dds)
resultsNames(ddstotal)
#[1] "Intercept"                 "species_AD4_vs_A2"         "species_D5_vs_A2"         
#[4] "species_PS7_vs_A2"         "species_TM1_vs_A2"         "condition_Salt_vs_Control"
#[7] "speciesAD4.conditionSalt"  "speciesD5.conditionSalt"   "speciesPS7.conditionSalt" 
#[10] "speciesTM1.conditionSalt"

###the treatment effect for each genotype
#The condition effect for A2 (the main effect for reference genotype)
results.resA2 = results(ddstotal, contrast=c("condition","Salt","Control"))
head(results.resA2)
print( summary(results.resA2,alpha=.05) ) #up-3324, down-3598
write.table(results.resA2,file="DEG.A2.txt",sep="\t")

#for the non-reference geneotype, adding the main effect and the interaction term
#The conidtion effect for D5 (the extra condition effect in D5 compared to A2)
results.resD5 = results(ddstotal, contrast = list(c("condition_Salt_vs_Control","speciesD5.conditionSalt")))
print(summary(results.resD5,alpha=.05)) #up-6063, down-6058
write.table(results.resD5,file="DEG.D5.txt",sep="\t")

#he condition effect for TM1 (the extra condition effect in TM1 compared to A2)
results.resTM1 = results(ddstotal, contrast = list(c("condition_Salt_vs_Control","speciesTM1.conditionSalt")))
print(summary(results.resTM1,alpha=.05)) #up-5864, down-5630
write.table(results.resTM1,file="DEG.TM1.txt",sep="\t")

#The condition effect for PS7 (the extra condition effect in PS7 compared to A2)
results.resPS7 = results(ddstotal,contrast = list(c("condition_Salt_vs_Control","speciesPS7.conditionSalt")))
print(summary(results.resPS7,alpha=.05)) #up-349, down-396
write.table(results.resPS7,file="DEG.PS7.txt",sep="\t")

#The condition effect for AD4 (the extra condition effect in AD4 compared to A2)
results.resAD4 = results(ddstotal, contrast = list(c("condition_Salt_vs_Control","speciesAD4.conditionSalt")))
print(summary(results.resAD4,alpha=.05,na.rm=TRUE)) #up-5240, down-4221
write.table(results.resAD4,file="DEG.AD4.txt",sep="\t")
#######################################################

###under control condition, the difference between D5,TM1,PS7, AD4 and A2: species
results.resD5A2.ck <- results(ddstotal, contrast = c(0,0,1,0,0,0,0,0,0,0))
head(results.resD5A2.ck)
print(summary(results.resD5A2.ck,alpha=.05,na.rm=TRUE)) #up-10328, down-9690
write.table(results.resD5A2.ck, file="DEG.D5vsA2ck.txt",sep="\t")

results.resTM1A2.ck <- results(ddstotal, contrast = c(0,0,0,0,1,0,0,0,0,0))
print(summary(results.resTM1A2.ck,alpha=.05,na.rm=TRUE)) #up-8199, down-7116
write.table(results.resTM1A2.ck, file="DEG.TM1vsA2ck.txt",sep="\t")

results.resPS7A2.ck <- results(ddstotal, contrast = c(0,0,0,1,0,0,0,0,0,0))
print(summary(results.resPS7A2.ck,alpha=.05,na.rm=TRUE)) #up-8092, down-7622
write.table(results.resPS7A2.ck, file="DEG.PS7vsA2ck.txt",sep="\t")

results.resAD4A2.ck <- results(ddstotal, contrast = c(0,1,0,0,0,0,0,0,0,0))
print(summary(results.resAD4A2.ck,alpha=.05,na.rm=TRUE)) #up-8221, down-7816
write.table(results.resAD4A2.ck, file="DEG.AD4vsA2ck.txt",sep="\t")

###under control condition, the difference between TM1,PS7, AD4 and D5
results.resTM1D5.ck <- results(ddstotal, contrast = c(0,0,-1,0,1,0,0,0,0,0))
print(summary(results.resTM1D5.ck,alpha=.05,na.rm=TRUE)) #up-9478, down-8923
write.table(results.resTM1D5.ck, file="DEG.TM1vsD5ck.txt",sep="\t")

results.resPS7D5.ck <- results(ddstotal, contrast = c(0,0,-1,1,0,0,0,0,0,0))
print(summary(results.resPS7D5.ck,alpha=.05,na.rm=TRUE)) #up-9042, down-9289
write.table(results.resPS7D5.ck, file="DEG.PS7vsD5ck.txt",sep="\t")

results.resAD4D5.ck <- results(ddstotal, contrast = c(0,1,-1,0,0,0,0,0,0,0))
print(summary(results.resAD4D5.ck,alpha=.05,na.rm=TRUE)) #up-9395, down-9580
write.table(results.resAD4D5.ck, file="DEG.AD4vsD5ck.txt",sep="\t")

###under salt treatment, the difference between D5, TM1, PS7, AD4 and A2
results.resD5A2.salt <- results(ddstotal, contrast = c(0,0,1,0,0,0,0,1,0,0))
print(summary(results.resD5A2.salt,alpha=.05,na.rm=TRUE)) #up-10760, down-9789
write.table(results.resD5A2.salt, file="DEG.D5vsA2salt.txt",sep="\t")

results.resTM1A2.salt <- results(ddstotal, contrast = c(0,0,0,0,1,0,0,0,0,1))
print(summary(results.resTM1A2.salt,alpha=.05,na.rm=TRUE)) #up-10148, down-8845
write.table(results.resTM1A2.salt, file="DEG.TM1vsA2salt.txt",sep="\t")

results.resPS7A2.salt <- results(ddstotal, contrast = c(0,0,0,1,0,0,0,0,1,0))
print(summary(results.resPS7A2.salt,alpha=.05,na.rm=TRUE)) #up-8883,down-8601
write.table(results.resPS7A2.salt, file="DEG.PS7vsA2salt.txt",sep="\t")

results.resAD4A2.salt <- results(ddstotal, contrast = c(0,1,0,0,0,0,1,0,0,0))
print(summary(results.resAD4A2.salt,alpha=.05,na.rm=TRUE)) #up-10153, down-9301
write.table(results.resAD4A2.salt, file="DEG.AD4vsA2salt.txt",sep="\t")

###under salt treatment, the difference between TM1, PS7, AD4 and D5
results.resTM1D5.salt <- results(ddstotal, contrast = c(0,0,-1,0,1,0,0,-1,0,1))
print(summary(results.resTM1D5.salt,alpha=.05,na.rm=TRUE)) #up-10363, down-9774
write.table(results.resTM1D5.salt, file="DEG.TM1vsD5salt.txt",sep="\t")

results.resPS7D5.salt <- results(ddstotal, contrast = c(0,0,-1,1,0,0,0,-1,1,0))
print(summary(results.resPS7D5.salt,alpha=.05,na.rm=TRUE)) #up-10036,down-11099
write.table(results.resPS7D5.salt, file="DEG.PS7vsD5salt.txt",sep="\t")

results.resAD4D5.salt <- results(ddstotal, contrast = c(0,1,-1,0,0,0,1,-1,0,0))
print(summary(results.resAD4D5.salt,alpha=.05,na.rm=TRUE)) #up-10660, down-9856
write.table(results.resAD4D5.salt, file="DEG.AD4vsD5salt.txt",sep="\t")

###interaction effect
#the interaction term for condition effect in D5, TM1, PS7, AD4 vs A2 
results.resD5A2int <- results(ddstotal,name="speciesD5.conditionSalt") #the difference in treatment effect btwn D5 and the reference genotype A2
results.resTM1A2int <- results(ddstotal, name="speciesTM1.conditionSalt")
results.resPS7A2int <- results(ddstotal, name="speciesPS7.conditionSalt")
resAD4A2int <- results(ddstotal, name="speciesAD4.conditionSalt")
#interaction term for condition effect in TM1, PS7, AD4 vs D5
results.resTM1D5int <- results(ddstotal, contrast = list("speciesTM1.conditionSalt","speciesD5.conditionSalt"))
results.resPS7D5int <- results(ddstotal, contrast = list("speciesPS7.conditionSalt","speciesD5.conditionSalt"))
results.resAD4D5int <- results(ddstotal, contrast = list("speciesAD4.conditionSalt","speciesD5.conditionSalt"))

###write output tables
#DEGs in response to salt treatment in each species
library("tidyverse")
DEG.A2.up <- results.resA2 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange > 0)
write.csv(DEG.A2.up, file="DEG_A2up.csv",row.names = F)

DEG.A2.down <- results.resA2 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange < 0)
write.csv(DEG.A2.down, file="DEG_A2down.csv",row.names = F)

DEG.D5.up <- results.resD5 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange > 0)
write.csv(DEG.D5.up, file="DEG_D5up.csv",row.names = F)

DEG.D5.down <- results.resD5 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange < 0)
write.csv(DEG.D5.down, file="DEG_D5down.csv",row.names = F)

DEG.TM1.up <- results.resTM1 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange > 0)
write.csv(DEG.TM1.up, file="DEG_TM1up.csv",row.names = F)

DEG.TM1.down <- results.resTM1 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange < 0)
write.csv(DEG.TM1.down, file="DEG_TM1down.csv",row.names = F)

DEG.PS7.up <- results.resPS7 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange > 0)
write.csv(DEG.PS7.up, file="DEG_PS7up.csv",row.names = F)

DEG.PS7.down <- results.resPS7 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange < 0)
write.csv(DEG.PS7.down, file="DEG_PS7down.csv",row.names = F)

DEG.AD4.up <- results.resAD4 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange > 0)
write.csv(DEG.AD4.up, file="DEG_AD4up.csv",row.names = F)

DEG.AD4.down <- results.resAD4 %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05 & log2FoldChange < 0)
write.csv(DEG.AD4.down, file="DEG_AD4down.csv",row.names = F)


###GO analysis of DE genes with topGO
library(topGO)
load('D5annotation.Rdata')
load('cottonGOenrich.RData')

#set up ID and GO term mapping
geneNames <- names(geneID2GO)

#make objects with DE results
deg.res = ls(pattern = 'results')
for(j in 1:length(deg.res)){
  deg.set <- get(deg.res[j])
  assign(paste(deg.res[j],".up",sep=""),deg.set[!is.na(deg.set$padj) & deg.set$log2FoldChange > 0,])
  assign(paste(deg.res[j],".down",sep=""),deg.set[!is.na(deg.set$padj) & deg.set$log2FoldChange < 0,])	
}

GOresults<-data.frame()

#set up table with different types of DGE terms
deg.full.res <- c(paste(deg.res,".up",sep=""),paste(deg.res,".down",sep=""))

#perform GO analysis over all types of DGE
for(m in 1:length(deg.full.res)){
  deg.set <- get(deg.full.res[m])
  sig.list <- rownames(deg.set[!is.na(deg.set$padj) & deg.set$padj < 0.05,])
  geneList <- factor(as.integer(geneNames %in% sig.list))
  if(length(levels(geneList)) != 2) next
  names(geneList) <- geneNames
  pdf(file=paste("GoTerms",deg.full.res[m],".pdf", sep=""))
  remove(enrich)
  for(on in c("MF","BP","CC"))
  {
    GOdata <- new("topGOdata", ontology = on, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    results.table <- GenTable(GOdata, result, topNodes = length(result@score))
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
    results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
    # label ontology type
    results.table$ontology<-on
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
    keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
    if(exists("enrich")) {enrich<- rbind(enrich, keep)}
    if(!exists("enrich")) {enrich<- keep}
    # draw figure for GO terms pval<=0.05 before FDR correction
    if(is.na(sigNo<-length(keep$ontology))){next}
    showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
    mtext(on, line=-1)
  }
  dev.off()
  #add enriched terms to whole GO enrichment results table
  if(dim(enrich)[1]>0)
  {
    enrichDE <- enrich
    enrichDE$DE <- deg.full.res[m]
    GOresults<-rbind(GOresults,enrichDE)
  }
}
write.table(GOresults, file="GOresults_deg_full.txt", sep="\t", row.names=FALSE)

###LRT tests for the effect of species, treatment and the interaction terms
#LRT to test effect of main terms
design(ddstotal) <- ~species + condition
ddsLRT.species <- DESeq(ddstotal, test="LRT", reduced = ~condition)
ddsLRT.condition <- DESeq(ddstotal, test="LRT", reduced = ~species)

DEG.species <- results(ddsLRT.species)
DEG.condition <- results(ddsLRT.condition)
DEG.condition.up.sig <- DEG.condition[!is.na(DEG.condition$log2FoldChange) & DEG.condition$log2FoldChange > 0 & !is.na(DEG.condition$padj) & DEG.condition$padj < 0.05,]
DEG.condition.down <- DEG.condition[!is.na(DEG.condition$log2FoldChange) & DEG.condition$log2FoldChange < 0,]

###the interactions term
dds1 <- DESeqDataSetFromMatrix(countData = counts.total[,-8], colData = coltotal29, design = ~ species + condition + species:condition)
LRT.dds <- DESeq(dds1,test="LRT", reduced = ~species + condition) 
DEG.interaction <- results(LRT.dds)

DEG.species.list <- DEG.species %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05) #30310
write.csv(DEG.species.list, file = "DEG_species.csv", row.names = F)

DEG.condition.list <- DEG.condition %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05) #13826
write.csv(DEG.condition.list, file = "DEG_condition.csv", row.names = F)

DEG.interaction.list <- DEG.interaction %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05) #15197
write.csv(DEG.interaction.list, file = "DEG_interaction.csv", row.names = F)

###GO enrichment
remove(enrich)
GOresults <- data.frame()

deg.red.res <- c("DEG.species","DEG.condition.up","DEG.condition.down","DEG.interaction")
for(m in 1:length(deg.red.res)){
  deg.set <- get(deg.red.res[m])
  sig.list <- rownames(deg.set[!is.na(deg.set$padj) & deg.set$padj < 0.05,])
  geneList <- factor(as.integer(geneNames %in% sig.list))
  if(length(levels(geneList)) != 2) next
  names(geneList) <- geneNames
  pdf(file=paste("GoTerms",deg.red.res[m],".pdf", sep=""))
  remove(enrich)
  for(on in c("MF","BP","CC"))
  {
    GOdata <- new("topGOdata", ontology = on, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    results.table <- GenTable(GOdata, result, topNodes = length(result@score))
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
    results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
    # label ontology type
    results.table$ontology<-on
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
    keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
    if(exists("enrich")) {enrich<- rbind(enrich, keep)}
    if(!exists("enrich")) {enrich<- keep}
    # draw figure for GO terms pval<=0.05 before FDR correction
    if(is.na(sigNo<-length(keep$ontology))){next}
    showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
    mtext(on, line=-1)
  }
  dev.off()
  #add enriched terms to whole GO enrichment results table
  if(dim(enrich)[1]>0)
  {
    enrichDE <- enrich
    enrichDE$DE <- deg.red.res[m]
    GOresults<-rbind(GOresults,enrichDE)
  }
}
write.table(GOresults, file="GOresults_deg_reduced.txt", sep="\t", row.names=FALSE)

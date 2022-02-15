### import counts data
counttable <- read.table("salt.counts",header=TRUE, sep="\t"); dim(counttable)
head(counttable);
names(counttable)
names(counttable) <- gsub("sort.bam","sort.T.bam",names(counttable))
names(counttable) <- gsub(".sort|.bam|_5.22|_5.23|_6.11","",names(counttable))
dim(saltcount <- counttable[,-c(2:6)])
rownames(saltcount) <- saltcount[,1]

# get subset data for each species
A2count <- saltcount[,grep("A2",names(saltcount))];dim(A2count)
D5count <- saltcount[,grep("D5",names(saltcount))];dim(D5count)
TM1count <- saltcount[,grep("TM1",names(saltcount))];dim(TM1count)
AD4count <- saltcount[,grep("AD4",names(saltcount))];dim(AD4count)
names(A2count)
# combine data from "web" in A2
A2ck1 <- A2count[,1:5]+A2count[,6:10]
A2ck3 <- A2count[,16:20]+A2count[,21:25]
A2salt2 <- A2count[,31:35]+A2count[,36:40]
A2salt3 <- A2count[,41:45]+A2count[,46:50]
A2count_c <- cbind(A2ck1,A2count[,grep("CK2",names(A2count))],A2ck3,A2count[,grep("Salt1",names(A2count))],A2salt2,A2salt3);dim(A2count_c)
counts <- cbind(A2count_c,D5count,TM1count,AD4count);dim(counts)

# pull out the Total counts values
counts.T <- counts[,grep(".T$",names(counts))];dim(counts.T)
head(counts.T)
colSums(counts.T)
save(counts.T, file = "counts.T.RData")

#pull out total counts of A2 and D5
diploid.T <- counts[,grep("A2.*T$|D5.*T$",names(counts))];dim(diploid.T)
#A and D counts of polyploids
TM1.AD <- counts[,grep("TM1.*A$|TM1.*D$",names(counts))];dim(TM1.AD)
#PS7.AD <- counts[,grep("PS7.*A$|PS7.*D$",names(counts))];dim(PS7.AD)
AD4.AD <- counts[,grep("AD4.*A$|AD4.*D$",names(counts))];dim(AD4.AD)

#get the count table for analysis, ADs subgenome
counts.sub <- cbind(diploid.T, TM1.AD, AD4.AD)
names(counts.sub) <- gsub("[.]T","",names(counts.sub))
save(counts.sub, file = "counts.sub.RData")

#At+Dt for allopolyploids
TM1.Total <- TM1.AD[,c(seq(1,11,by=2))]+TM1.AD[,c(seq(2,12,by=2))]
AD4.Total <- AD4.AD[,c(seq(1,11,by=2))]+AD4.AD[,c(seq(2,12,by=2))]
counts.total <- cbind(diploid.T, TM1.Total, AD4.Total)
names(counts.total) <- gsub("[.]T|[.]A","",names(counts.total))
save(counts.total, file = "counts.total.RData")

## loading packages
library(DESeq2)
library(ggplot2)
library(scales)
library(ape)
library(pheatmap)
library(RColorBrewer)

# make a coldata
species <- c(rep("A2",6),rep("D5",6),rep("AD1",6),rep("AD4",6))
condition <- rep(c(rep("Control",3),rep("Salt",3)),4)
rep <- rep(1:3,8)
coltotal <- data.frame(row.names=names(counts.total),species,condition,rep)
coltotal$sample <- paste(coltotal$species,coltotal$condition,sep = ".")
all(colnames(counts.total)==rownames(coltotal))
coltotal$condition <- relevel(coltotal$condition,ref = "Control")
write.table(coltotal, file="coldatatotal.txt",row.names = T)

# Build DESeq data
dds <- DESeqDataSetFromMatrix(countData = counts.total, colData = coltotal, design = ~ species + condition + species:condition)

# rlog normalization
rld <- rlog(dds, blind = FALSE)
count.rld <- as.data.frame(assay(rld))

mytheme <- theme_bw()+
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

pdf(file="PCA_rld.pdf")
ggplot(aes(PC1, PC2, color=species, shape=condition),data=dat) + geom_point(size=5) +xlab(proportion[1]) + ylab(proportion[2])+mytheme
dev.off()

#heatmap
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
head(rld_cor)
pdf(file="correlation_rld.pdf")
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
#[1] "Intercept"                 "species_AD1_vs_A2"         "species_AD4_vs_A2"        
#[4] "species_D5_vs_A2"          "condition_Salt_vs_Control" "speciesAD1.conditionSalt" 
#[7] "speciesAD4.conditionSalt"  "speciesD5.conditionSalt"

###the treatment effect for each genotype
#The condition effect for A2 (the main effect for reference genotype)
results.resA2 = results(ddstotal, contrast=c("condition","Salt","Control"))
head(results.resA2)
print( summary(results.resA2,alpha=.05) ) #up-2875, down-3220
write.table(results.resA2,file="DEG.A2.txt",sep="\t")

#for the non-reference geneotype, adding the main effect and the interaction term
#The conidtion effect for D5 (the extra condition effect in D5 compared to A2)
results.resD5 = results(ddstotal, contrast = list(c("condition_Salt_vs_Control","speciesD5.conditionSalt")))
print(summary(results.resD5,alpha=.05)) #up-5830, down-5714
write.table(results.resD5,file="DEG.D5.txt",sep="\t")

#he condition effect for AD1 (the extra condition effect in TM1 compared to A2)
results.resTM1 = results(ddstotal, contrast = list(c("condition_Salt_vs_Control","speciesAD1.conditionSalt")))
print(summary(results.resTM1,alpha=.05)) #up-5534, down-5156
write.table(results.resTM1,file="DEG.TM1.txt",sep="\t")

#The condition effect for AD4 (the extra condition effect in AD4 compared to A2)
results.resAD4 = results(ddstotal, contrast = list(c("condition_Salt_vs_Control","speciesAD4.conditionSalt")))
print(summary(results.resAD4,alpha=.05,na.rm=TRUE)) #up-4883, down-3923
write.table(results.resAD4,file="DEG.AD4.txt",sep="\t")

### write output tables
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

### GO analysis of DE genes with topGO
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

###the interactions term
dds1 <- DESeqDataSetFromMatrix(countData = counts.total[,-8], colData = coltotal29, design = ~ species + condition + species:condition)
LRT.dds <- DESeq(dds1,test="LRT", reduced = ~species + condition) 
DEG.interaction <- results(LRT.dds)

DEG.species.list <- DEG.species %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05) #29054
write.csv(DEG.species.list, file = "DEG_species.csv", row.names = F)

DEG.condition.list <- DEG.condition %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05) #13378
write.csv(DEG.condition.list, file = "DEG_condition.csv", row.names = F)

DEG.interaction.list <- DEG.interaction %>%
  data.frame() %>%
  rownames_to_column(var="Geneid") %>% 
  as_tibble()%>%
  filter(padj < 0.05) #12979
write.csv(DEG.interaction.list, file = "DEG_interaction.csv", row.names = F)

###GO enrichment
remove(enrich)
GOresults <- data.frame()

deg.red.res <- c("DEG.species","DEG.condition","DEG.interaction")
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
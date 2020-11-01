######################
##Load the packages###
######################
setwd("D://inISU/SaltRNA/WGCNA_results/")
library(WGCNA)
library(RColorBrewer)
library(flashClust)
library(stringr)
library(ggplot2)
library(BiocParallel)
library(DESeq2)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

################################
##Read in the expression data###
################################

counttable <- read.csv("counts.total.txt",sep="")
colinfo <- read.csv("coldatatotal.txt",sep="")

# Make each row corresponds to a gene and column to a sample or auxiliary information.
saltdat = as.data.frame(t(counttable));
write.table(t(saltdat),file="W-saltcount.raw.txt", sep="\t")

pdf("W-raw.boxplot.pdf")
boxplot(log2(t(saltdat) ), las=2)
dev.off()

# Normalization using DESeq2 rlog
dds_log <- DESeqDataSetFromMatrix( countData = counttable, colData = colinfo, design = ~ sample)
# rlog transformation, note that with defalt blind=TRUE, design is actually ~1
rld <- rlog(dds_log)
head(assay(rld))
pdf("s1.rld_pca.pdf")
plotPCA(rld, intgroup = c("species", "condition"))
dev.off()

expr<-as.data.frame(assay(rld) )
names(expr)<-names(counttable)
expr<-expr[,-which(names(expr) %in% "D5_CK2")] #D5_20 is a sample from polyploid, remove it
write.table(expr,"W-count.rlog.txt", sep="\t")

###Read in the physiological traits data
traitData<-read.table("saltTraits.txt", header=TRUE, sep="\t")
rownames(traitData) = traitData[, 1]
datTraits <- traitData[,-1] #for condition, 0 represents control,and 1 reps salt stress
# make sure the rownames in datTraits line up with those in expression data
datTraits <- datTraits[-8,]
table(rownames(datTraits)==names(expr))

# construct expr
datExprT<-t(as.data.frame(expr))

# sample network based on squared Euclidean distance
# note that we transpose the data
A=adjacency(expr,type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)

# Designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-5 # often -2.5

# the color vector indicates outlyingness (red)
outlierColor1=ifelse(Z.k<thresholdZ.k,"red","black")

thresholdZ.k= -2.5
outlierColor2=ifelse(Z.k<thresholdZ.k ,"red","black")

pdf("s1.sample_dendrogram_and_trait_heatmap.rlog.pdf")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation:
# where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits),"C",sep="")
datColors=data.frame(outlierC_5=outlierColor1,outlierC_2.5=outlierColor2, traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap")
dev.off()

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenes(datExprT, verbose = 3);
gsg$allOK
# Excluding 1019 genes from the calculation due to too many missing samples or zero variance.

# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExprT)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExprT)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExprT = datExprT[gsg$goodSamples, gsg$goodGenes]
}

save(datExprT, datTraits, file = "R-01-dataInput.RData")

########################
##Choosing softPower####
########################

#set strings as factors to false (as in WGCNA documentation)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
lnames = load(file = "R-01-dataInput.RData")
lnames
nGenes<-ncol(datExprT) #36486L

# Choose a set of soft-thresholding powers
type<-"signed"
powers = c(c(1:10), seq(from = 12, to=40, by=2))

sft.signed = pickSoftThreshold(datExprT, powerVector = powers, networkType=type,verbose = 5)
#Basically, the default power is 6 for unsigned network, 16 for signed network; if the observed power choice is smaller than default, I will choose the observed power, otherwise use default.
sft.signed$powerEstimate #24 

# Plot the results:
pdf("s2.ChooseSoftThresholdPower_16.pdf")
par(mfrow = c(1,2));
cex1 = 0.7;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft.signed$fitIndices[,1], -sign(sft.signed$fitIndices[,3])*sft.signed$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));

text(sft.signed$fitIndices[,1], -sign(sft.signed$fitIndices[,3])*sft.signed$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft.signed$fitIndices[,1], sft.signed$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft.signed$fitIndices[,1], sft.signed$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

###########################
###Network Construction####
###########################

power = 24
#cor <- WGCNA::cor #avoid the conflict between WGCNA function with the R base stat, can be reset to cor<-stats::cor after finished
#This step is quite slow and takes a lot of memory, so I actually ran it on a server; otherwise it may split into multiple blocks but difference were found between one single block and multiple blocks
net = blockwiseModules(datExprT, 
                       checkMissingData = TRUE,
                       blocks = NULL, randomSeed = 12345,
                       maxBlockSize = nGenes,
                       power = power, 
                       networkType = "signed",TOMType = "signed", corType = "pearson",
                       mergeCutHeight = 0.25,
                       deepSplit = 2,
                       
                       minModuleSize = min(30, ncol(datExprT)/2 ),
                       
                       pamStage = TRUE, pamRespectsDendro = TRUE,
                       
                       reassignThreshold = 0,
                       numericLabels = TRUE,verbose = 3,
                       
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Allspecies_24_power_signed_merge_0.25_TOM" 
)
assign("saltnet24", net)
save("saltnet24", file = "R-02-buildNetwork24.RData")
net=saltnet24
table(net$colors)

###plot the modules and dendrogram
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
pdf("S3.module_colors_dendrogram.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

##save the results
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
# get eigengenes
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "R-03-allspecies_networkConstruction-auto.RData")

MEs_col = MEs 
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
#######################################
###General network topology analysis###
#######################################
load("R-01-dataInput.RData")
load("R-02-buildNetwork24.RData")
net = saltnet24
# Displaying module heatmap and the eigengene
samples <- c(rep("A2_CK",3),rep("A2_Salt",3),rep("D5_CK",2),rep("D5_Salt",3),rep("TM1_CK",3),rep("TM1_Salt",3),rep("PS7_CK",3),rep("PS7_Salt",3),rep("AD4_CK",3),rep("AD4_Salt",3))
ss <- as.factor(samples)

Nmodules = dim(net$MEs)[2]
MEs<-net$MEs

subDat   <-  apply(datExprT,2,as.numeric)
plots <- list()  # new empty list
pdf("s4.eigengene_expression_heatmap.pdf")
for(me in 0:(Nmodules-1)){
  
  which.module=paste("ME",me,sep="")
  module.color=labels2colors(me)
  #heatmap
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
  plotMat(t(scale(subDat[,net$colors==me ]) ),
          nrgcols=30,rlabels=T,rcols=module.color,
          main=paste(which.module, module.color, sep=": "), cex.main=2)
  #barplot
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(MEs[,which.module], col=module.color, main="", cex.main=2,
          ylab="eigengene expression",xlab="treatment", names.arg=as.character(ss))
  
}
dev.off()

pdf("s5.eigengene_expression_association_with_traits.pdf")
for(me in 0:(Nmodules-1)){
  which.module=paste("ME",me,sep="")
  module.color=labels2colors(me)
  #line, anova
  df<-data.frame(ME=MEs[,which.module], ss, module = which.module )
  fit<-aov(ME~ss,df)
  dfc<-summarySE(df, measurevar="ME", groupvars=c("ss", "module"))
  dfc$species <- factor(gsub("_CK|_Salt","",dfc$ss),levels=c("A2","D5","TM1","PS7","AD4"))
  dfc$condition <- gsub("A2_|D5_|TM1_|PS7_|AD4_","",dfc$ss)
  plots[[me+1]]<- ggplot(dfc, aes(x=species, y=ME, fill = condition)) +
    geom_bar(stat="identity",position=position_dodge(), color="black", size=0.3) +
    geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.3,position=position_dodge(0.9)) +
    ggtitle(paste(which.module," ",module.color, sep="") )+
    theme_bw() +
    theme(plot.title=element_text( size=11),legend.position = "none")
  
  for(page in 1:ceiling(Nmodules/9))
  {
    if(Nmodules>(9*page))
    {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
    else
    {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
  }
}
dev.off()

##########################################################
###Relate modules to phenotype and functional gene sets###
##########################################################

#load("R-02-buildNetwork.RData")
# eigengene~sample, anova
pval<-apply(MEs,2,function(x){round(anova(aov(x~ss) )$"Pr(>F)"[1],4)})
pval<-as.data.frame(pval)
pval$symbol<-ifelse(pval$pval<0.05,"*"," ")
pval$numeric<-as.numeric(substring(rownames(pval),3) )
pval<-pval[order(pval$numeric),]
pval$symbol[1]<-" "  # ME0 always meaningless
pval

## or run linear model to test for population, treatment, and interaction effects for module expression for each module
#borrowed from Mead et al., 2019 Mol. Eol.,  https://github.com/alaynamead/valley_oak_water_stress
species <- c(rep("A2",6),rep("D5",6),rep("TM1",6),rep("PS7",6),rep("AD4",6))
condition <- rep(c(rep("Control",3),rep("Salt",3)),5)
rep <- rep(1:3,10)
coltotal <- data.frame(row.names=names(datExprT),species,condition,rep)
col29 <- coltotal[-8,]
row.names(col29) <- row.names(net$MEs)
col29$sample <- row.names(net$MEs)

net$MEs <- net$MEs[,order(colnames(net$MEs))]
head(net$MEs)
lm.results <- list()
for(m in 1:ncol(net$MEs)){
  lm.results[[m]] <- list(anova(lm(net$MEs[,m] ~ col29$species*col29$condition)))
}
names(lm.results) <- colnames(MEs)
save(lm.results, file = 'R-00-lm_results_module_expression_sitextreatment_power24.Rdata')

# make GxE plots with average expression level for each site
col29$species <- factor(col29$species, levels = c('A2', 'D5', 'TM1', 'PS7', 'AD4'))
col29$species
# species colors
cols <- c("#FDD692","#67D5B5","#EE7785","#C89EC4", "#000262")
# function for GxE plots
gxe_plot <- function(expr.df, mod, species = c('A2', 'D5', 'TM1', 'PS7', 'AD4'), cols=c("#FDD692","#67D5B5","#EE7785","#C89EC4", "#000262"), main = mod, legendsize = 0.8, legend = T, lty = rep(1,6), lwd = 6, ...){
  ex <- expr.df[,mod]
  names(ex) <- rownames(expr.df)
  ex.species <- as.data.frame(matrix(nrow = 5, ncol = 6))
  rownames(ex.species) <- species
  colnames(ex.species) <- c('mean.c', 'mean.d', 'median.c', 'median.d', 'sd.c', 'sd.d')
  
  for(n in 1:length(species)){
    # get list of control samples for each species
    samples.c <- col29$sample[col29$species == species[n] & col29$condition == 'Control']
    ex.species[n,1] <- mean(as.numeric(ex[as.character(samples.c)]))
    ex.species[n,3] <- median(as.numeric(ex[as.character(samples.c)]))
    ex.species[n,5] <- sd(as.numeric(ex[as.character(samples.c)]))
    
    # get list of salt samples for each species
    samples.s <- col29$sample[col29$species == species[n] & col29$condition == 'Salt'] 
    ex.species[n,2] <- mean(as.numeric(ex[as.character(samples.s)]))
    ex.species[n,4] <- median(as.numeric(ex[as.character(samples.s)]))
    ex.species[n,6] <- sd(as.numeric(ex[as.character(samples.s)]))
  }
  
  min = min(c(ex.species$mean.c, ex.species$mean.d) - 0.05*abs(max(c(ex.species$mean.c, ex.species$mean.d) - min(c(ex.species$mean.c, ex.species$mean.d)))))
  max = max(c(ex.species$mean.c, ex.species$mean.d) + 0.05*abs(max(c(ex.species$mean.c, ex.species$mean.d) - min(c(ex.species$mean.c, ex.species$mean.d)))))
  
  plot(c(ex.species$mean.c[1], ex.species$mean.d[1]), type = 'l', ylim = c(min,max), col = cols[1], xlab = '', ylab = 'Module Expression', main = main, xaxt="n", lty = 1, lwd = lwd, ...)
  axis(1, labels = c('Control', 'Salt'), at = c(1,2))
  lines(c(ex.species$mean.c[2], ex.species$mean.d[2]), col = cols[2], lwd = lwd, lty = 2)
  lines(c(ex.species$mean.c[3], ex.species$mean.d[3]), col = cols[3], lwd = lwd, lty = 3)
  lines(c(ex.species$mean.c[4], ex.species$mean.d[4]), col = cols[4], lwd = lwd, lty = 4)
  lines(c(ex.species$mean.c[5], ex.species$mean.d[5]), col = cols[5], lwd = lwd, lty = 5)
  if(legend == T){
    legend(x = 'top', fill = cols, legend = species, horiz = T, cex = legendsize)
  }
}
par(mfrow = c(1,1), mar = c(5,4,4,1)+0.1, las = 1)
###
plot.window.orig <- plot.window
# plot single module:
gxe_plot(expr.df = net$MEs, mod = 'ME1', species = levels(col29$species), lty = c(1:6),  main = 'ME1\n(description)')
#plot GxE plots for all modules
for(n in 1:ncol(net$MEs)){
  mod <- colnames(net$MEs)[n]
  png(filename = paste('module_GxE_', mod, '.png', sep = ''), res = 300, width = 7, height = 6, units = 'in')
  par(mfrow = c(1,1), mar = c(5,4,4,1)+0.1, las = 1)
  gxe_plot(expr.df = net$MEs, mod = mod, species = levels(col29$species))
  dev.off()
}

## get list of modules with significant species, condition, and interaction effects
mod.s <- vector()
mod.c <- vector()
mod.i <- vector()
for(n in 1:length(lm.results)){
  
  # species
  if(lm.results[[n]][[1]][1,5] <= 0.05){
    mod.s <- append(mod.s, names(lm.results)[n])
  }
  # condition
  if(lm.results[[n]][[1]][2,5] <= 0.05){
    mod.c <- append(mod.c, names(lm.results)[n])
  }
  # interaction
  if(lm.results[[n]][[1]][3,5] <= 0.05){
    mod.i <- append(mod.i, names(lm.results)[n])
  }
}
mod.s 
mod.c 
mod.i  
# remove grey, not a real module
mod.c <- mod.c[which(mod.c != 'grey')]
# plot all significant modules
length(unique(c(mod.s, mod.i, mod.c)))
mods <- unique(c(mod.s, mod.i, mod.c))

## plot heatmap of average expression differences between control/treatment for each population
# data
col29$interaction <- interaction(col29$condition, col29$species)
# make df with avg module expression by site and treatment
me.avg <- as.data.frame(matrix(nrow = ncol(net$MEs), ncol = 10))
colnames(me.avg) <- levels(col29$interaction)
rownames(me.avg) <- colnames(net$MEs)
for(mod in 1:ncol(net$MEs)){
  
  for(group in 1:length(levels(col29$interaction))){
    g <- levels(col29$interaction)[group]
    me.avg[mod, group] <- mean(net$MEs[which(col29$interaction == g), mod])
  }
}
# use only treatment and interaction modules
me.avg <- me.avg[unique(c(mod.i, mod.c)),]
# looking at treatment differences
colnames(me.avg)
me.avg.diff <- me.avg[,c(2,4,6,8,10)] - me.avg[,c(1,3,5,7,9)]
colnames(me.avg.diff) <- gsub('Salt.', '', colnames(me.avg.diff))
cols <- colorRamp2(breaks = as.vector(quantile(unlist(me.avg.diff))), colors = c('darkslateblue', 'lightcyan2', 'white', 'mistyrose','firebrick3'))
# to include functions
#rownames(me.avg.diff) <- c("blue\n(chloroplast)\nPxT *", "darkgreen\n(ribosome)\nT***, PxT**", "greenyellow\n(DNA replication)\nPxT*", "midnightblue\n¨D\nPxT*", "pink\n(kinase activity)\nT**, PxT*", "yellow\n(oxidoreductase activity)\nP*, PxT*", "black\n(response to stress)\nT***", "brown\n¨D\nT*", "grey60\n(protein folding)\nP*, T***")
# to add population colors
cols.top <- as.list(c("#FDD692","#67D5B5","#EE7785","#C89EC4", "#000262"))
names(cols.top) <- colnames(me.avg.diff)
top <- HeatmapAnnotation(Population = colnames(me.avg.diff), col = list(Population = c('MC' =  '#ffe55c', 'FT' = '#ff8b4e', 'FH' = '#e14d66', 'CV' = '#9c2e7f', 'PL' = '#5f0092', 'RD' = '#000262')), annotation_legend_param = list(title='\nPopulation'))
#png(file = 'plots/module_expression_avg_change_heatmap_with_functions.png', res = 300, height = 5.5, width = 7.5, units = 'in')
Heatmap(me.avg.diff,
        cluster_columns = T,
        cluster_rows = T,
        col = cols,
        row_names_side = 'left',
        #column_names_side = 'top',
        heatmap_legend_param = list(title = 'Average\nExpression\nDifference'),
        #top_annotation = top,
        column_title_gp = gpar(fontsize = 14))
#dev.off()


###Eigengene adjacency heatmap
pdf("s6.Eigengene_adjacency_heatmap_1.pdf",width = 12,height = 12)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), cex.lab = 0.8,
                      plotDendrograms = T,
                      xLabelsAngle = 90)
dev.off()


## add the species info to traits data
datTraits_spe <- datTraits
datTraits_spe$species_A2 <- c(rep(1,6),rep(0,23))
datTraits_spe$species_D5 <- c(rep(0,6),rep(1,5),rep(0,18))
datTraits_spe$species_TM1 <- c(rep(0,11),rep(1,6),rep(0,12))
datTraits_spe$species_PS7 <- c(rep(0,17),rep(1,6),rep(0,6))
datTraits_spe$species_AD4 <- c(rep(0,23),rep(1,6))

### Relate eigengenes to external traits or sample conditions
MET=orderMEs(cbind(MEs,datTraits_spe))
#Visualization of the eigengene network representing the relationships among the modules and sample traits. The top panel shows a hierarchical clustering dendrogram of the eigengenes based on the dissimilarity diss(q_1,q_2)=1-cor(E^{(q_1)},E^{(q_2)}). The bottom panel shows the shows the eigengene adjacency A_{q1,q2}=0.5+0.5 cor(E^{(q_1)},E^{(q_2)}).
plotEigengeneNetworks(MET,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)
# plot eigengene network for only the significant modules
module.sig<-rownames(pval[pval$symbol=="*",])
MET.sig<-MET[,module.sig]
plotEigengeneNetworks(MET.sig,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)
# plot it again with color names
names(MET.sig)<-paste("ME",labels2colors(as.numeric(substring(names(MET.sig),3) ) ),sep="")
plotEigengeneNetworks(MET.sig,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)

### graphical representation for correlation with modules
#get correlation and p-vals
moduleTraitCor = cor(MEs, datTraits_spe, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples=nrow(datTraits_spe));
MEcolors<-paste("ME",labels2colors(as.numeric(gsub("ME","",names(MEs)) ) ), sep="")
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1), mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
# Table of module-trait correlations and p-values. Each cell reports the correlation (and p-value) resulting from  correlating module eigengenes (rows) to traits (columns). The table is color-coded by correlation according to the color legend.
pdf("s7.ModuleTraitAssociation.pdf", width=16, height=16)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(moduleTraitCor), yLabels = MEcolors, ySymbols = names(MEs), colorLabels = TRUE, colors = blueWhiteRed(50), textMatrix = as.matrix(textMatrix), setStdMargins = FALSE, cex.text = 0.5,zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

# plot it for only significant modules
where.sig<-sort(match(module.sig, rownames(moduleTraitCor)) )
moduleTraitCor.sig <- moduleTraitCor[where.sig,]
textMatrix.sig <- textMatrix[where.sig,]
labeledHeatmap(Matrix = moduleTraitCor.sig, xLabels = colnames(moduleTraitCor), yLabels = MEcolors[where.sig], ySymbols = rownames(moduleTraitCor.sig), colorLabels = TRUE, colors = blueWhiteRed(50), textMatrix = as.matrix(textMatrix.sig), setStdMargins = FALSE, cex.text = 0.7,zlim = c(-1,1), main = paste("Module-trait relationships: sig only"))

# For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.
# calculate the module membership values
# (aka. module eigengene based connectivity kME):
MM=signedKME(datExprT, MEs)
rownames(MM)<-names(datExprT)

# define a gene significance variable. We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait.
GS=cor(datExprT,datTraits_spe,use="p")
# This translates the numeric values into colors
GS.Color=numbers2colors(GS,signed=T)

# Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules. As an example, we look at the brown module that a high correlation with body weight. We plot a scatterplot of Gene Significance vs. Module Membership in select modules...
colorOfColumn=substring(names(MM),4)
pdf("s8.GSvsMM.pdf")
par(mfrow = c(3,1))
for (module in sigModule) {
  column = match(module,colorOfColumn)
  restModule=net$colors==module
  for (trait in colnames(GS)){
    verboseScatterplot(MM[restModule,column], GS[restModule,trait], xlab=paste("Module Membership of ME",module),ylab=paste("GS.",trait),    main=paste("kME.",module,"vs. GS"), col=labels2colors(module))
  }
}
dev.off()


#######################
###make the TOM plot###
#######################
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^16
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
pdf("s9.TOM_plot.pdf")
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")  
dev.off()


##########################
### Modules Annotation ###
##########################

####### write gene with corresponding module assignment and annotation
aa<-load('D5annotation.Rdata')
aa  
aa<-annotation221[,c("transcript","tair10.defline", "pfam","panther","kog", "kegg.ec", "kegg.orthology", "go", "tair10", "tair10.symbol")]
dim(aa<-aa[grep("[.]1$",aa$transcript),])  #37505
aa$gene = gsub("[.]1$","", aa$transcript)
me <- as.data.frame(net$colors)
names(me) <- "ME"
me$gene <- rownames(me)
rownames(me) <- NULL
dim(me<-merge(me, aa, all.x=TRUE, by="gene"))
write.table(me, file="s5.module&annotation.txt", row.names=FALSE,sep="\t")


######## Functional enrichment analysis
# add Gorai id to module group
MM$ID<-colnames(datExprT)
#for(me in rownames(pval))
#{
#  me<-paste("k",me,sep="")
#  write.table(MM[,c("ID",me)], file=paste("GSEA_MM_rnk",me,".rnk",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
#}
# Run GseaPreranked analysis
#MM<-MM[,order(names(MM))]
#write.table(MM,file="s5.moduleMembership.txt",sep="\t",row.names=FALSE,quote=FALSE)
#save(saltnet16, datTraits_spe, pval, GS, MM, file = "R-05-all12_GS&MM.RData")

me0<-me
names(MM)[1]<-"gene"
dim(me0<-merge(me0,MM,by="gene",all.x=TRUE, all.y=TRUE) )
# me0[me0$kME9>0.9 & me0$ME==9, c("gene","tair10.defline")]
write.table(me0,file="s5.moduleMembership&annotation.txt",sep="\t",row.names=FALSE)


#######Annotation with topGO
library(topGO)
load('D5annotation.Rdata')
BiocManager::install("GSEABase")
load('cottonGOenrich.RData')
# "geneID2GO"       "go4profile"  "goraiUniverse"   "gsc"             "level1_BP_terms" "level1_MF_terms"  "level2_BP_terms" "level2_MF_terms" "level3_BP_terms" "level3_MF_terms"   "level4_BP_terms" "level4_MF_terms" "level5_BP_terms" "level5_MF_terms"
# all included in individual network
universe<-colnames(datExprT)
# for each module containing genes of interest
GOresults<-data.frame()
for(module in 0:(Nmodules-1))
{
  genes<-universe[net$colors==module]
  geneList <- factor(as.integer(universe %in% genes))
  names(geneList) <- universe
  
  pdf(file=paste("ME",module,".pdf", sep=""))
  
  # topGO analysis
  remove(enrich)
  for(on in c("MF","BP","CC"))
  {
    print(on)
    # Make topGO object
    GOdata <- new("topGOdata", ontology = on, allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    # fisher test
    result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    results.table <- GenTable(GOdata, result, topNodes = length(result@score))
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
    results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
    # label ontology type
    results.table$ontology<-on
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
    keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
    if(exists("enrich")) enrich<- rbind(enrich, keep)
    if(!exists("enrich")) enrich<- keep
    
    # draw figure for GO terms pval<=0.05 before FDR correction
    if(is.na(sigNo<-length(keep$ontology))){next}
    showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
    mtext(on, line=-1)
  }
  dev.off()
  if(dim(enrich)[1]>0)
  {
    enrichME<-enrich
    enrichME$ME=module
    GOresults<-rbind(GOresults,enrichME)   }
}
write.table(GOresults, file="s8.GOresults.txt", sep="\t", row.names=FALSE)
####################################
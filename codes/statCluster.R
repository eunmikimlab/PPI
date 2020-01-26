#--------------------------------------------------------------------------------------------------
# Purpose: Draw MTB HeatMap for PPI Motif frequency
# Date	 : Sep 24, 2018
# Author : Jong Cheol Jeong (JongCheol.Jeong@uky.edu) 
#--------------------------------------------------------------------------------------------------
rm(list=ls()); graphics.off(); closeAllConnections();

#------------------------------------------------------------------------------
# Environment Setup
#------------------------------------------------------------------------------
reqPackages <- c('shiny', 'maps', 'mapproj', 'RMySQL','data.table','ggplot2','reshape2', 'RColorBrewer', 'randomcoloR','DT', 'plyr', 'BBmisc', 'devtools')
newPackages <- reqPackages[!(reqPackages %in% installed.packages()[,'Package'])]
if (length(newPackages) > 0) {
  install.packages(newPackages, repos='http://cran.r-project.org', dependencies=TRUE);
}

#// for Bioconductor
reqPackages <- c("ComplexHeatmap")
newPackages <- reqPackages[!(reqPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(newPackages, dependencies=TRUE)  
}

#// using source codes
reqPackages <- c("circlize")
newPackages <- reqPackages[!(reqPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) {
  install.packages('http://cran.rstudio.com/src/contrib/circlize_0.4.1.tar.gz', repo=NULL, dependencies=TRUE, type='source');
}

library('shiny')
library('maps')
library('mapproj')
library('RMySQL')
library('data.table')
library('ggplot2')
library('reshape2')
library('ComplexHeatmap')
library('circlize')
library('RColorBrewer')
library('randomcoloR')
library('DT')
library('plyr')
library('BBmisc')
library('dendextend')

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
options(expressions=500000)
dbName <- 'EC'
len_motif <- 1
inFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/heatmap_%s_seqFreq_lenMotif_%d.csv', toupper(dbName), len_motif)
outFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/heatmap_%s_seqFreq_lenMotif_%d.pdf', toupper(dbName), len_motif)
pcaFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/pca_%s_seqFreq_lenMotif_%d.pdf', toupper(dbName), len_motif)

inData <- read.csv(inFile, header=TRUE, sep=",")

inData$NAME <- aaply(as.vector(unlist(inData$CHAIN)), 1, function(x) {tmp = strsplit(x, '[.]')[[1]][1]})
inData$GROUP1  <- aaply(as.vector(unlist(inData$CATEGORY)), 1, function(x) {tmp = strsplit(x, '[.]')[[1]][1]})
inData$GROUP2  <- aaply(as.vector(unlist(inData$CATEGORY)), 1, function(x) {tmp = strsplit(x, '[.]')[[1]]; paste(tmp[1], tmp[2], sep='.')})

Names <- as.vector(unique(unlist(inData$GROUP2)))

#define heatmap (genotype) colors 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colVector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
muCol     <- c('ghostwhite', colVector[1:(length(Names)-1)])
muColors  <- structure(muCol, names=Names)


#Draw HeatMap
dataMtx <- inData[, 1:20]
dataMtx <- dataMtx[order(inData$GROUP2), ]
inData$INDEX <- 1:dim(dataMtx)[1]
#rownames(dataMtx) <- inData$INDEX

#dataMtx <- as.data.frame(lapply(dataMtx, function(x) (x-min(x))/(max(x)-min(x))))
#dataMtx <- scale(t(dataMtx), center=TRUE, scale=TRUE);
#rownames(normMtx) <- inData$NAME

inData$NEWNAME <- paste(inData$GROUP2, inData$INDEX, sep = '_')
rownames(dataMtx) <- inData$NEWNAME

num_cluster <- length(unique(inData$GROUP2))
dend <- hclust(dist(dataMtx))
#dend <- color_branches(dend, k = num_cluster)
Clusters <- cutree(dend, k = num_cluster)
cIndex <- as.vector(Clusters)
cValue <- aaply(names(Clusters), 1, function(x) {tmp <- strsplit(strsplit(x, '_')[[1]][1], '[.]')[[1]]; as.numeric(paste(tmp, collapse = ''))})
cor.test(cIndex, cValue, method=c('spearman'))


### PCA
groups <- as.vector(unlist(inData$GROUP1))

#-- run  PCA and plot the significance of PCs
dataPCA <- prcomp(dataMtx, center=TRUE, scale=TRUE)
plot(dataPCA, type='l')
summary(dataPCA)

library(ggbiplot)

pdf(pcaFile, width=12, height=10)
par(cex=0.4, mar=c(3,1,1,15));
g <- ggbiplot(dataPCA, choices=c(1,2), obs.scale=1, var.scale=1, 
              groups=groups, ellipse=TRUE, var.axes=FALSE, 
              #labels=groups, 
              circle=TRUE)
g <- g + scale_color_discrete(name = 'Date group')
g <- g + theme(legend.direction = 'vertical', legend.position = 'right')
print(g)
dev.off()
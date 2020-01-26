#--------------------------------------------------------------------------------------------------
# Purpose: Draw MTB HeatMap for PPI Motif frequency
# Date	 : Sep 24, 2018
# Author : Jong Cheol Jeong (JongCheol.Jeong@uky.edu) 
#--------------------------------------------------------------------------------------------------
rm(list=ls()); graphics.off(); closeAllConnections();

#------------------------------------------------------------------------------
# Environment Setup
#------------------------------------------------------------------------------
reqPackages <- c('shiny', 'maps', 'mapproj', 'RMySQL','data.table','ggplot2','reshape2', 'RColorBrewer', 'randomcoloR','DT', 'plyr', 'BBmisc', 'devtools', 'rafalib')
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
library(ape)

options(expressions=500000)
#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
dbName <- 'EC'
selBase <- 'seqBase'     # seqBase=sequence-based interface residues, strcBase=structure-based interface residues
selIntf <- 'NoIntf'         # target to search motif: 'Intf'= interface residue, 'NoIntf'=noninterface residues
len_motif <- 1

  
cmt <- sprintf('[%s] %s with %s and len_motif %s is processing.. ', dbName, selBase, selIntf, len_motif)
write(cmt, '')

if (selBase == 'seqBase') {
  inFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/heatmap_%s_seqFreq_lenMotif_%d.csv',selBase, selIntf, toupper(dbName), len_motif)
  outFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/heatmap_%s_seqFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  inFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/heatmap_%s_strcFreq_lenMotif_%d.csv',selBase, selIntf, toupper(dbName), len_motif)
  outFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/heatmap_%s_strcFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
}

# read data 
#----------------------------------------------------------
inData <- read.csv(inFile, header=TRUE, sep=",")

#-- convert all NAs to 0
inData[is.na(inData)] <- 0

#-- add Names for analysis into extra columns
inData$NAME <- aaply(as.vector(unlist(inData$CHAIN)), 1, function(x) {tmp = strsplit(x, '[.]')[[1]][1]})
inData$GROUP1  <- aaply(as.vector(unlist(inData$CATEGORY)), 1, function(x) {tmp = strsplit(x, '[.]')[[1]][1]})
inData$GROUP2  <- aaply(as.vector(unlist(inData$CATEGORY)), 1, function(x) {tmp = strsplit(x, '[.]')[[1]]; paste(tmp[1], tmp[2], sep='.')})

# Define row names
#----------------------------------------------------------
col_data <- dim(inData)[2] - 8 
dataMtx <- inData[, 1:col_data]
dataMtx <- dataMtx[order(inData$GROUP2), ]
inData$INDEX <- 1:dim(dataMtx)[1]

# GROUP 1 is the highest EC hierarchy
# GROUP is the second highest EC hierarchy
Names <- inData$GROUP2
inData$NEWNAME <- paste(Names, inData$INDEX, sep = '_')
rownames(dataMtx) <- inData$NEWNAME

# Box plot of data distribution
#----------------------------------------------------------
st_dataMtx <- stack(dataMtx)

if (selBase == 'seqBase') {
  outBox <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/boxData_%s_seqFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  outBox <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/boxData_%s_strcFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
}
pdf(outBox, width=12, height=10)
boxplot(st_dataMtx$values ~ st_dataMtx$ind)
dev.off()


# Build Hierarchical cluster
# https://www.rdocumentation.org/packages/rafalib/versions/1.0.0
#----------------------------------------------------------
hc = hclust(dist(dataMtx))
num_cluster <- length(unique(Names))

if (selBase == 'seqBase') {
  outHcl <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/hclust_%s_seqFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  outHcl <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/hclust_%s_strcFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
}
pdf(outHcl, width=15, height=200)
#myplclust(hc, labels=Names, lab.col=as.fumeric(Names), cex=0.2, hang=-1)
plot(as.phylo(hc), cex=0.3)
dev.off()

# Get Hierarchical Information 
#----------------------------------------------------------
if (selBase == 'seqBase') {
  infoHcl <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/hclustInfo_%s_seqFreq_lenMotif_%d.csv', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  infoHcl <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/hclustInfo_%s_strcFreq_lenMotif_%d.csv', selBase, selIntf, toupper(dbName), len_motif)  
}
listHc <- as.data.frame(cutree(hc, k=100))
colnames(listHc) <- c('CLUSTID')
listHc$SMPNAME <- rownames(listHc)

write.csv(listHc, infoHcl, row.names=FALSE)

# Get cluster information
#----------------------------------------------------------
Clusters <- cutree(hc, k=num_cluster)
cIndex <- as.vector(Clusters)
#cValue <- aaply(names(Clusters), 1, function(x) {tmp <- strsplit(strsplit(x, '_')[[1]][1], '[.]')[[1]]; as.numeric(paste(tmp, collapse = ''))})
cValue <- aaply(names(Clusters), 1, function(x) {strsplit(x, '_')[[1]][1]})

#-- store frequency table
if (selBase == 'seqBase') {
  outConf <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/Confusion_%s_seqFreq_lenMotif_%d.csv', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  outConf <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/Confusion_%s_strcFreq_lenMotif_%d.csv', selBase, selIntf, toupper(dbName), len_motif)  
}

freMtx <- table(cIndex, cValue)
dfMtx <- data.frame(unclass(freMtx))
colnames(dfMtx) <- sort(unique(Names))
write.csv(dfMtx, outConf)

# ANOVA test
# http://statland.org/R/R/R1way.htm
#----------------------------------------------------------
# Make each column is a clust
tMtx <- t(dfMtx)
rownames(tMtx) <- rownames(freMtx)
dftMtx <- data.frame(unlist(tMtx))
colnames(dftMtx) <- rownames(freMtx)
#colnames(dftMtx) <- paste(rownames(freMtx), 'cl', sep='_')
stMtx <- stack(dftMtx)
summary(stMtx)
if (selBase == 'seqBase') {
  outClust <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/boxClust_%s_seqFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  outClust <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/boxClust_%s_strcFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
}

pdf(outClust, width=12, height=10)
boxplot(stMtx$values ~ stMtx$ind)
dev.off()

if (selBase == 'seqBase') {
  outAnova <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/Anova_%s_seqFreq_lenMotif_%d.txt', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  outAnova <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/Anova_%s_strcFreq_lenMotif_%d.txt', selBase, selIntf, toupper(dbName), len_motif)  
}

outS <- oneway.test(stMtx$values ~ stMtx$ind, var.equal=TRUE)
capture.output(outS, file=outAnova)
#aov(values ~ ind)
#cor.test(cIndex, cValue, method=c('spearman'))


# PCA
#----------------------------------------------------------
#-- define how many groups needed to be clustered 
groups <- as.vector(unlist(inData$GROUP2))

#-- run  PCA and plot the significance of PCs
#dataPCA <- prcomp(dataMtx, center=TRUE, scale=TRUE)
dataPCA <- prcomp(dataMtx, center=TRUE)
if (selBase == 'seqBase') {
  outEffPca <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/pcaEffect_%s_seqFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  outEffPca <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/pcaEffect_%s_strcFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
}
pdf(outEffPca, width=12, height=10)
plot(dataPCA, type='l')
dev.off()

if (selBase == 'seqBase') {
  outSpca <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/PCA_%s_seqFreq_lenMotif_%d.txt', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  outSpca <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/PCA_%s_strcFreq_lenMotif_%d.txt', selBase, selIntf, toupper(dbName), len_motif)  
}

Spca <- summary(dataPCA)
capture.output(Spca, file=outSpca)

library(ggbiplot)
if (selBase == 'seqBase') {
  pcaFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/pcaClust_%s_seqFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  pcaFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/pcaClust_%s_strcFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
}

pdf(pcaFile, width=12, height=10)
par(cex=0.4, mar=c(3,1,1,15));
g <- ggbiplot(dataPCA, choices=c(1,2), obs.scale=1, var.scale=1, 
              groups=groups, ellipse=TRUE, var.axes=FALSE, 
              #labels=groups, 
              circle=TRUE)
g <- g + scale_color_discrete(name = 'Clusters')
g <- g + theme(legend.direction = 'vertical', legend.position = 'right')
print(g)
dev.off()
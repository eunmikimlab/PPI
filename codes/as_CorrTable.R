#--------------------------------------------------------------------------------------------------
# Purpose: Draw MTB HeatMap for PPI Motif frequency
# Date	 : Sep 24, 2018
# Author : Jong Cheol Jeong (JongCheol.Jeong@uky.edu) 
#--------------------------------------------------------------------------------------------------
rm(list=ls()); graphics.off(); closeAllConnections();

#------------------------------------------------------------------------------
# Environment Setup
#------------------------------------------------------------------------------
reqPackages <- c('shiny', 'maps', 'mapproj', 'RMySQL','data.table','ggplot2','reshape2', 'RColorBrewer', 'randomcoloR','DT', 'plyr', 'BBmisc', 'devtools', 'rafalib', 'corrplot')
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
library('ape')
library('corrplot')

options(expressions=500000)
#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
dbName <- 'SCOP'
selBase <- 'strcBase'     # seqBase=sequence-based interface residues, strcBase=structure-based interface residues
selIntf <- 'Intf'         # target to search motif: 'Intf'= interface residue, 'NoIntf'=noninterface residues
len_motif <- 3
#height <- 0.2
num_clust <- 100

if (len_motif == 1) {
  h <- 10
  w <- 10
} else {
  h <- 30
  w <- 30
}

  
cmt <- sprintf('[%s] %s with %s and len_motif %s is processing.. ', dbName, selBase, selIntf, len_motif)
write(cmt, '')

if (selBase == 'seqBase') {
  inFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/heatmap_%s_seqFreq_lenMotif_%d.csv',selBase, selIntf, toupper(dbName), len_motif)
  outFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/corTable_%s_seqFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  inFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/heatmap_%s_strcFreq_lenMotif_%d.csv',selBase, selIntf, toupper(dbName), len_motif)
  outFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/corTable_%s_strcFreq_lenMotif_%d.pdf', selBase, selIntf, toupper(dbName), len_motif)  
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

M <- cor(dataMtx, use="complete.obs")
M[is.na(M)] <- 0

pdf(outFile, width=w, height=h)
corrplot(M, addrect = 2, order = "hclust", method = "ellipse")
dev.off()


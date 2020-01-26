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


# Get Hierarchical Information 
#----------------------------------------------------------
hc = hclust(dist(dataMtx))

if (selBase == 'seqBase') {
  infoHcl <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/hclustInfo_%s_seqFreq_lenMotif_%d.csv', selBase, selIntf, toupper(dbName), len_motif)  
} else if (selBase == 'strcBase') {
  infoHcl <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/hclustInfo_%s_strcFreq_lenMotif_%d.csv', selBase, selIntf, toupper(dbName), len_motif)  
}
#listHc <- as.data.frame(cutree(hc, h=height))
listHc <- as.data.frame(cutree(hc, k=num_clust))
colnames(listHc) <- c('CLUSTID')
listHc$SMPNAME <- rownames(listHc)

num_cluster = max(unique(listHc$CLUSTID))
#cmt <- sprintf('With H=%f total %d clusters are identified', height, num_cluster)
#print(cmt)

# Identifying UniProt ID 
# The results will be upload to https://www.uniprot.org/uploadlists/
# IDs may not be completly matched
#----------------------------------------------------------
if (selBase == 'seqBase') {
  outFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/UniProt_%s_seqFreq_lenMotif_%d_CL%d.csv', selBase, selIntf, toupper(dbName), len_motif, num_clust)  
} else if (selBase == 'strcBase') {
  outFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/%s/%s/UniProt_%s_strcFreq_lenMotif_%d_CL%d.csv', selBase, selIntf, toupper(dbName), len_motif, num_clust)  
}

infoMtx <- merge(listHc, inData, by.x='SMPNAME', by.y='NEWNAME', all=TRUE)[, c('SMPNAME', 'PDB', 'NAME')]
tmp <- aaply(infoMtx$NAME, .margins=1, .fun=function(x){strsplit(x, '_')[[1]][3]})
infoMtx$CHAIN <- as.vector(tmp)

if (dbName == 'EC') {
  dbParse <- '/Users/jjeong/local/project_dev/ppi/dbs/pdb_chain_enzyme.csv'
  dfParse <- read.csv(dbParse, header=TRUE, comment.char = "#")
  mergeMtx <- merge(infoMtx, dfParse, by.x=c('PDB', 'CHAIN'), by.y=c('PDB', 'CHAIN'))
  fullMtx <- merge(mergeMtx, listHc, by.x=c('SMPNAME'), by.y=c('SMPNAME'))
  
} else if (dbName == 'SCOP') {
  dbParse = '/Users/jjeong/local/project_dev/ppi/dbs/pdb_chain_scop_uniprot.csv'
  dfParse <- read.csv(dbParse, header=TRUE, comment.char = "#")
  mergeMtx <- merge(infoMtx, dfParse, by.x=c('PDB', 'CHAIN'), by.y=c('PDB', 'CHAIN'))
  fullMtx <- merge(mergeMtx, listHc, by.x=c('SMPNAME'), by.y=c('SMPNAME'))
  
}

write.csv(fullMtx, file=outFile, row.names = FALSE)


#--------------------------------------------------------------------------------------------------
# Purpose: Draw MTB HeatMap for PPI Motif frequency
# Date	 : Sep 24, 2018
# Author : Jong Cheol Jeong (JongCheol.Jeong@uky.edu) 
#--------------------------------------------------------------------------------------------------
rm(list=ls()); graphics.off(); closeAllConnections();

#------------------------------------------------------------------------------
# Environment Setup
#------------------------------------------------------------------------------
reqPackages <- c('shiny', 'maps', 'mapproj', 'RMySQL','data.table','ggplot2','reshape2', 'RColorBrewer', 'randomcoloR','DT', 'plyr', 'BBmisc')
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
outFile <- sprintf('/Users/jjeong/local/project_dev/ppi/outputs/heatmap/heatmap_%s_seqFreq_lenMotif_%d_norm.pdf', toupper(dbName), len_motif)
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

dend = hclust(dist(dataMtx))
dend = color_branches(dend, k = length(unique(inData$GROUP2)))

pdf(outFile, width=15, height=200)

Heatmap(dataMtx, 
                name = 'Motif Frequency',
                column_names_side = 'top',
                row_names_side = 'left',
                row_dend_side = 'right',
                cluster_rows = dend,
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                show_row_dend = TRUE,
                row_names_gp = gpar(fontsize = 3),
                row_dend_width = unit(50, "mm")
              )
dev.off()


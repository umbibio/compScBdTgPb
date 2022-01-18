library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(parallel)
library(openxlsx)

source('./util_funcs.R')


num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
intra.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.intra.expr.csv"

## IDs
prod.desc  <- read.xlsx('../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx')


getExpr <- function(in.file, TGGT1_ME49){
  file.counts <- read.csv(in.file)
  genes <- file.counts$X
  ind <- which(genes %in% TGGT1_ME49$TGME49)
  file.counts <- file.counts[ind, ]
  genes <- genes[ind]
  genes <- TGGT1_ME49$TGGT1[match(genes, TGGT1_ME49$TGME49)]
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  return(expr)
}

intra.counts <- getExpr(intra.file.csv, TGGT1_ME49)


## individual Seurat objects
S.O.intra <- CreateSeuratObject(counts = intra.counts)
S.O.intra$orig.ident <- 'intra'
set.seed(100)
S.O.intra <- subset(x = S.O.intra, downsample = 6000)

feats <- c("nFeature_RNA","nCount_RNA")
VlnPlot(S.O.intra,features = feats, pt.size = 0.1,ncol = 4) + NoLegend()
FeatureScatter(S.O.intra, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)
S.O.intra <- subset(x = S.O.intra, subset = nFeature_RNA > 80 & nFeature_RNA < 1400 )
dim(S.O.intra)



## Individually process the data and transfer labels

## Boothroyed data
S.O.tg <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.RData')
S.O.intra <- prep_S.O(S.O.intra)
anchors <- FindTransferAnchors(reference = S.O.tg, query = S.O.intra, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = S.O.tg@meta.data$phase,dims = 1:30)
predictions$phase <- predictions$predicted.id
S.O.intra <- AddMetaData(object = S.O.intra, metadata = predictions)


## For updating the APP
saveRDS(S.O.intra, '../scExpressionProfiler/S_O_toxo_MJ_labels.rds')



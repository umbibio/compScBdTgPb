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
library(sctransform)
library(glassoFast)
library(igraph)
library(ggraph)
library(graphlayouts)
library(Matrix)


source('util_funcs.R')

## Estimate empirical covariance matrix from expression data.
## Network smoothing is applied to each expression matrix with parameters alpha = 0.3 and max.iter = 3
## All data are down sampled to include 800 cells/cluster max and 2000 most variable features.
## The network in k-nn extracted from the Seurat object.


### b. divergens
input.dir.bdiv <- "../Input/scRNAseqBdiv/"
bdiv.count.file <- "bdiv.expr.csv"
## Reading b. divergens data
bdiv.count <- read.csv(paste(input.dir.bdiv, bdiv.count.file, sep = ''))
genes <- bdiv.count$X
bd.expr <- bdiv.count[,-1]
rownames(bd.expr) <- genes
# Set initial Seurat clusters
S.O.bd <- CreateSeuratObject(counts = bd.expr, min.cells = 10, min.features = 100)
S.O.bd <- subset(S.O.bd, subset = nFeature_RNA > 200 & nFeature_RNA < 1200 )
S.O.bd <- prep_S.O(S.O.bd, res = 0.1, var.features = T, down.sample = T, smooth.data = T, network = T)
Ft.bd  <- S.O.bd@assays$smooth@data

sds <- apply(Ft.bd, 1, sd)
rm.ind <- which(sds == 0)

if(length(rm.ind) > 0){
  Ft.bd <- Ft.bd[-rm.ind,]
}

Ft.bd.scale <- scale(t(Ft.bd))
S.bd <- cov(Ft.bd.scale)

saveRDS(S.bd, '../Input/compScBdTgPb/glasso/S_bd_smooth_2000.RData')

### T. gondii
input.dir.tg <- '../Input/scRNAseqToxoAtlas/kz/'
rh384.expr.file <- 'rh384_expression_filtered.csv'
## Reading T. gondii data
rh384.expr <- read.csv(paste(input.dir.tg, rh384.expr.file, sep = ''))

processCounts <- function(expr){
  cols <- expr[,1]
  expr <- t(expr[,-1])
  colnames(expr) <- cols
  return(expr)
}

rh384.expr <- processCounts(rh384.expr)
# Set initial Seurat clusters
S.O.tg <- CreateSeuratObject(counts = rh384.expr, min.cells = 10, min.features = 100)
#VlnPlot(S.O.tg, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O.tg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.tg <- subset(S.O.tg, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
S.O.tg <- prep_S.O(S.O.tg, res = 0.1, var.features = T, down.sample = T, smooth.data = T, network = T)
Ft.tg <- S.O.tg@assays$smooth@data
Ft.tg.scale <- scale(t(Ft.tg))
S.tg <- cov(Ft.tg.scale)
saveRDS(S.tg, '../Input/compScBdTgPb/glasso/S_tg_smooth_2000.RData')

## Plasmodium
input.dir.10x <- "../Input/MalariaCellAtlas/Expression_Matrices/10X/pb10xIDC/"
plasmodium.10x.count.file <- 'pb10xIDC_counts.csv'
plasmodium.10x.count <- read.csv(paste(input.dir.10x, plasmodium.10x.count.file, sep = ''))
genes.10x <- plasmodium.10x.count$X
pb.10x.expr <- plasmodium.10x.count[,-1]
rownames(pb.10x.expr) <- genes.10x
S.O.pb <- CreateSeuratObject(counts = pb.10x.expr, min.cells = 10, min.features = 100)
VlnPlot(S.O.pb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.pb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.pb <- subset(S.O.pb, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 )
S.O.pb <- prep_S.O(S.O.pb, res = 0.1, var.features = T, down.sample = T, smooth.data = T, network = T)
Ft.pb <- S.O.pb@assays$smooth@data
Ft.pb.scale <- scale(t(Ft.pb))
S.pb <- cov(Ft.pb.scale)
saveRDS(S.pb, '../Input/compScBdTgPb/glasso/S_pb_smooth_2000.RData')




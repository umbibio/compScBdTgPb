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
### T. gondii
S.O.bd.gam <- readRDS('../Input/compScBdTgPb/RData/S.O.bd.gam.RData')
Ft.bd <- S.O.bd.gam@assays$smooth@data
sds <- apply(Ft.bd, 1, sd)
rm.ind <- which(sds == 0)

if(length(rm.ind) > 0){
  Ft.bd <- Ft.bd[-rm.ind,]
}

## Time-course reconstruction
bd_cell.cycle.genes.df <- readRDS('../Input/compScBdTgPb/RData/bd_cell_cycle_genes_df.RData')

bd_sc.tc.df <- bd_cell.cycle.genes.df %>% 
  transmute(GeneID = GeneID, Sample = Sample, t = t, rep = rep)

Ft.bd_long <- Ft.bd %>% data.frame() %>% mutate(GeneID = rownames(.)) %>% 
  pivot_longer(-GeneID, names_to = 'Sample', values_to = 'raw_expr') 

Ft.bd_long$GeneID <- gsub('-', '_', Ft.bd_long$GeneID)

Ft.bd_long <- left_join(Ft.bd_long, bd_sc.tc.df, by = c('GeneID', 'Sample'))

Ft.bd_long.fits <- Ft.bd_long %>% group_by(GeneID, t) %>% summarise(mean_trend = mean(raw_expr))
Ft.bd_long <- left_join(Ft.bd_long, Ft.bd_long.fits, by = c('GeneID', 't'))
Ft.bd_long <- Ft.bd_long %>% arrange(GeneID, t, rep)
Ft.bd_long <- Ft.bd_long %>% mutate(de_trend.exp = raw_expr - mean_trend)

gg <- unique(Ft.bd_long$GeneID)
par(mfrow=c(2,1))
for(i in 1:20){
  tmp <- Ft.bd_long %>% dplyr::filter(GeneID == gg[i])
  plot(tmp$t, tmp$raw_expr)
  points(tmp$t, tmp$mean_trend, type = 'l', col = 'red')
  plot(tmp$t, tmp$de_trend.exp)
  Sys.sleep(1)

}

Ft.bd_wide <- Ft.bd_long %>% dplyr::select(Sample, GeneID, de_trend.exp) %>%
  pivot_wider(everything(), names_from = Sample, values_from = de_trend.exp) 

genes <- gsub('_', '-', Ft.bd_wide$GeneID)
Ft.bd_wide <- Ft.bd_wide %>% dplyr::select(-GeneID)
rownames(Ft.bd_wide) <- genes

Ft.bd.scale <- scale(t(Ft.bd))
S.bd <- cov(Ft.bd.scale)
saveRDS(S.bd, '../Input/compScBdTgPb/glasso/S_bd_smooth_gam_de_trended.RData')



### T. gondii
S.O.tg.gam <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.gam.RData')
Ft.tg <- S.O.tg.gam@assays$smooth@data
sds <- apply(Ft.tg, 1, sd)
rm.ind <- which(sds == 0)

if(length(rm.ind) > 0){
  Ft.tg <- Ft.tg[-rm.ind,]
}

## Time-course reconstruction
tg_cell.cycle.genes.df <- readRDS('../Input/compScBdTgPb/RData/tg_cell_cycle_genes_df.RData')

tg_sc.tc.df <- tg_cell.cycle.genes.df %>% 
  transmute(GeneID = GeneID, Sample = Sample, t = t, rep = rep)

Ft.tg_long <- Ft.tg %>% data.frame() %>% mutate(GeneID = rownames(.)) %>% 
  pivot_longer(-GeneID, names_to = 'Sample', values_to = 'raw_expr') 

Ft.tg_long$GeneID <- gsub('-', '_', Ft.tg_long$GeneID)

Ft.tg_long <- left_join(Ft.tg_long, tg_sc.tc.df, by = c('GeneID', 'Sample'))

Ft.tg_long.fits <- Ft.tg_long %>% group_by(GeneID, t) %>% summarise(mean_trend = mean(raw_expr))
Ft.tg_long <- left_join(Ft.tg_long, Ft.tg_long.fits, by = c('GeneID', 't'))
Ft.tg_long <- Ft.tg_long %>% arrange(GeneID, t, rep)
Ft.tg_long <- Ft.tg_long %>% mutate(de_trend.exp = raw_expr - mean_trend)

# gg <- unique(Ft.tg_long$GeneID)
# par(mfrow=c(2,1))
# for(i in 1:20){
#   tmp <- Ft.tg_long %>% dplyr::filter(GeneID == gg[i])
#   plot(tmp$t, tmp$raw_expr)
#   points(tmp$t, tmp$mean_trend, type = 'l', col = 'red')
#   plot(tmp$t, tmp$de_trend.exp)
#   Sys.sleep(1)
#   
# }

Ft.tg_wide <- Ft.tg_long %>% dplyr::select(Sample, GeneID, de_trend.exp) %>%
  pivot_wider(everything(), names_from = Sample, values_from = de_trend.exp) 

genes <- gsub('_', '-', Ft.tg_wide$GeneID)
Ft.tg_wide <- Ft.tg_wide %>% dplyr::select(-GeneID)
rownames(Ft.tg_wide) <- genes

Ft.tg.scale <- scale(t(Ft.tg))
S.tg <- cov(Ft.tg.scale)
saveRDS(S.tg, '../Input/compScBdTgPb/glasso/S_tg_smooth_gam_de_trended.RData')

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




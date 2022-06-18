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
library(cowplot)
library(patchwork)
library(doParallel)


source('./util_funcs.R')

## Extracting the data for control over PCA directions
getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), phase = S.O@meta.data$phase, 
                          spp = S.O@meta.data$spp)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)



## IDs
prod.desc  <- read.xlsx('../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx')


## Boothroyd data from Pru streesed parasites.
count.file <- "../Input/scRNAseqToxoAtlas/kz/pru_expression_filtered.csv"
meta.file  <- "../Input/scRNAseqToxoAtlas/kz/pru_obs_filtered.csv"
meta.file2  <- "../Input/compScBdTgPb/boothroyd_sc_all_data/Data/pru/adata_allg2/obs.csv" ## python generated

pru.counts <- read.csv(count.file)
#pru.meta <- read.csv(meta.file)
pru.meta <- read.csv(meta.file2)
pru.meta$cell_cycle <- gsub('G1 b', 'G1.b', gsub('G1 a', 'G1.a', gsub("\"", "", pru.meta$cell_cycle)))
pru.meta$dpi <- gsub(" ", "", pru.meta$dpi)
pru.meta$renamed_cluster

samples <- pru.counts$index
genes <- colnames(pru.counts)
ind <- which(genes %in% TGGT1_ME49$TGME49)
pru.counts <- pru.counts[, ind]
genes <- genes[ind]
genes <- TGGT1_ME49$TGGT1[match(genes, TGGT1_ME49$TGME49)]
colnames(pru.counts) <- genes

d0.ind <- pru.meta$dpi == 'Day0'
d3.ind <- pru.meta$dpi == 'Day3'
d5.ind <- pru.meta$dpi == 'Day5'
d7.ind <- pru.meta$dpi == 'Day7'

d0.expr <- data.frame(t(pru.counts[samples %in% pru.meta$index[d0.ind], ]))
colnames(d0.expr) <- samples[samples %in% pru.meta$index[d0.ind]]
d0.meta <- pru.meta[d0.ind, ] %>% transmute(Sample = index, dpi = dpi, phase = cell_cycle, cluster = renamed_cluster)
rownames(d0.meta) <- d0.meta$Sample
S.O.d0 <- CreateSeuratObject(counts = d0.expr)
S.O.d0$orig.ident <- 'pru.d0'
S.O.d0 <- AddMetaData(S.O.d0, d0.meta)


d3.expr <- data.frame(t(pru.counts[samples %in% pru.meta$index[d3.ind], ]))
colnames(d3.expr) <- samples[samples %in% pru.meta$index[d3.ind]]
d3.meta <- pru.meta[d3.ind, ] %>% transmute(Sample = index, dpi = dpi, phase = cell_cycle,  cluster = renamed_cluster)
rownames(d3.meta) <- d3.meta$Sample
S.O.d3 <- CreateSeuratObject(counts = d3.expr)
S.O.d3$orig.ident <- 'pru.d3'
S.O.d3 <- AddMetaData(S.O.d3, d3.meta)

d5.expr <- data.frame(t(pru.counts[samples %in% pru.meta$index[d5.ind], ]))
colnames(d5.expr) <- samples[samples %in% pru.meta$index[d5.ind]]
d5.meta <- pru.meta[d5.ind, ] %>% transmute(Sample = index, dpi = dpi, phase = cell_cycle, cluster = renamed_cluster)
rownames(d5.meta) <- d5.meta$Sample
S.O.d5 <- CreateSeuratObject(counts = d5.expr)
S.O.d5$orig.ident <- 'pru.d5'
S.O.d5 <- AddMetaData(S.O.d5, d5.meta)


d7.expr <- data.frame(t(pru.counts[samples %in% pru.meta$index[d7.ind], ]))
colnames(d7.expr) <- samples[samples %in% pru.meta$index[d7.ind]]
d7.meta <- pru.meta[d7.ind, ] %>% transmute(Sample = index, dpi = dpi, phase = cell_cycle, cluster = renamed_cluster)
rownames(d7.meta) <- d7.meta$Sample
S.O.d7 <- CreateSeuratObject(counts = d7.expr)
S.O.d7$orig.ident <- 'pru.d7'
S.O.d7 <- AddMetaData(S.O.d7, d7.meta)

S.O.list <- list(d0 = S.O.d0, d3 = S.O.d3, d5 = S.O.d5, d7 = S.O.d7)

saveRDS(S.O.list, '../Input/compScBdTgPb/RData/S.O.tg.pru.boothroyd.list2.rds')

S.O.merge <- merge(S.O.list[[1]], S.O.list[2:4], add.cell.ids=c("pru.d0","pru.d3","pru.d5","pru.d7"))

#### Quick test (no ancholing)
S.O.merge <- prep_S.O(S.O.merge)
S.O.merge@meta.data$spp <- S.O.merge@meta.data$dpi

Idents(S.O.merge) <- 'cluster'
DimPlot(S.O.merge, reduction = 'pca')
DimPlot(S.O.merge, reduction = 'umap')

SAG1.ID <- 'TGGT1-233460'
BAG1.ID <- "TGGT1-259020"


p <- FeaturePlot(S.O.merge, features = c(SAG1.ID, BAG1.ID), split.by = 'spp')
plot(p)

S.Os.list <- SplitObject(S.O.merge, split.by = "spp")
S.Os.list <- lapply(X = S.Os.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = S.Os.list)
anchors <- FindIntegrationAnchors(object.list = S.Os.list, anchor.features = features)
S.Os.integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(S.Os.integrated) <- "integrated"
S.Os.integrated <- ScaleData(S.Os.integrated, verbose = FALSE)
S.Os.integrated <- RunPCA(S.Os.integrated, npcs = 30, verbose = FALSE)
S.Os.integrated <- RunUMAP(S.Os.integrated, reduction = "pca", dims = 1:30)
S.Os.integrated <- FindNeighbors(S.Os.integrated, reduction = "pca", dims = 1:30)
S.Os.integrated <- FindClusters(S.Os.integrated, resolution = 0.2)
Idents(S.Os.integrated) <- 'cluster'
DimPlot(S.Os.integrated, reduction = 'umap', label = T)

DefaultAssay(S.Os.integrated) <- "RNA"
p <- FeaturePlot(S.Os.integrated, features = c(SAG1.ID, BAG1.ID))
plot(p)


## Clusters 7 & 8 correspond to original cluster P6
## Other markers
my.genes <- c("TGGT1-208020",
              "TGGT1-240900",
              "TGGT1-280460",
              "TGGT1-306620",
              "TGGT1-214960",
              "TGGT1-203050",
              "TGGT1-309410",
              "TGGT1-264485",
              "TGGT1-269010",
              "TGGT1-215895",
              'TGGT1-200385',
              "TGGT1-259020",
              'TGGT1-233460')


gene.names <- c("AP2Ib-1",
                "AP2VI-2",
                "AP2VIIa-2",
                'AP2IX-9',
                "AP2X-8",
                "AP2VIIa-6",
                "AP2XI-1",
                "AP2IX-3",
                "AP2VIII-7", 
                "AP2IX-10",
                "BFD1",
                "BAG1",
                "SAG1")


#i <- 3
i <- 10
FeaturePlot(S.Os.integrated, features = my.genes[i])



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
library(rhdf5)


source('./util_funcs.R')
set.seed(100)


getS.O <- function(count.mat, spp,  min.feature = 100, max.feature = 1200){
  S.O <- CreateSeuratObject(counts = count.mat, min.cells = 10, min.features = 100)
  S.O[["percent.mt"]] <- PercentageFeatureSet(object = S.O, pattern = "^BmR1-mt")
  #VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #FeatureScatter(object = Bm, feature1 = "nCount_RNA", feature2 = "percent.mt") 
  #FeatureScatter(object = Bm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

  S.O <- subset(x = S.O, subset = nFeature_RNA > min.feature & nFeature_RNA < max.feature &  percent.mt < 3)
  #S.O  <- subset(S.O , subset = nFeature_RNA > 100 & nFeature_RNA < 1200)

  
  pheno <- data.frame(Sample = names(S.O$orig.ident))
  pheno$spp <- spp
  pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = '_')
  
  L <- list(pheno = pheno, S.O = S.O)
  
  return(L)
}

## Read B microti data
h5_data <- hdf5r::H5File$new("../Input/Bmicroti/filtered_feature_bc_matrix.h5", mode = 'r')

feature_matrix <- Matrix::sparseMatrix(
  i = h5_data[['matrix/indices']][],
  p = h5_data[['matrix/indptr']][],
  x = h5_data[['matrix/data']][],
  dimnames = list(
    h5_data[['matrix/features/name']][],
    h5_data[['matrix/barcodes']][]
  ),
  dims = h5_data[['matrix/shape']][],
  index1 = FALSE
)

expr.bm <- as.matrix(feature_matrix)

## Read B. Div
file.counts <- read.csv("../Input/scRNAseqBdivCS/bd0hr_count.csv")
genes <- file.counts$X
expr.bd <- file.counts[,-1]
rownames(expr.bd) <- genes


## Get orthologous genes
Bd.Bm.orth <- read.xlsx('../Input/orthologs/Bdiv_Bmic_orth.xlsx')
Bd.Bm.orth <- inner_join(inner_join(Bd.Bm.orth, 
                                    data.frame(Bd = rownames(expr.bd)), by = 'Bd'), 
                         data.frame(Bm = rownames(expr.bm)), by = 'Bm')

##Get Sexual Stage gene
Sextual.Stage <- read.xlsx('../Input/SextualStageMarkers/Babesia_Sexual_Stage_Markers_KZ.xlsx', sheet = 2)


keep.bd <- which(rownames(expr.bd) %in% Bd.Bm.orth$Bd)
expr.bd <- expr.bd[keep.bd, ]

keep.bm <- which(rownames(expr.bm) %in% Bd.Bm.orth$Bm)
expr.bm <- expr.bm[keep.bm, ]
#rownames(expr.bd) <- Bd.Bm.orth$Bm[match(rownames(expr.bd), Bd.Bm.orth$Bd)]
rownames(expr.bm) <- Bd.Bm.orth$Bd[match(rownames(expr.bm), Bd.Bm.orth$Bm)]


S.O.bm <- getS.O(expr.bm, 'Bm', 50, 250)

S.O.bdiv <- getS.O(expr.bd, 'Bdiv', 100, 1200)


S.O.integrated <- processeMergedS.O(list(S.O.bdiv, S.O.bm))
S.O.integrated@meta.data$spp <- factor(S.O.integrated@meta.data$spp,
                                       levels = c('Bdiv', 'Bm'))

#prod.desc <- read.csv('../Input/BdivCellCycle/ProductDescription/Bm_Prod_desc.csv')
prod.desc <- read.csv('../Input/BdivCellCycle/ProductDescription/BDvi_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

#Idents(all.samples.integrated) <- "phase.cond"
p1 <- DimPlot(S.O.integrated, reduction = "pca", dims = c(1,2),
              #group.by = "cells", 
              split.by = 'spp',
              pt.size = 1,
              shape.by='spp',
              label = TRUE, label.size = 6) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p1)

ggsave(filename="../Output/scClockFigs/b_microti_b_divergence_pca.pdf",
       plot=p1,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600)

## Individual plots
DefaultAssay(S.O.integrated) <- "RNA"
Idents(S.O.integrated) <- "spp"
Idents(S.O.integrated) <- "phase.cond"

my.gene <- "Bdiv-038820"
sextual.marker <- gsub(" ", "", gsub('_', '-', Sextual.Stage$Bdiv))
p2 <- FeaturePlot(object = S.O.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = sextual.marker,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)

S.O.list <- SplitObject(S.O.integrated, split.by = "spp")

p3 <- FeaturePlot(object = S.O.list[[2]], 
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = sextual.marker[c(1,3)], blend = TRUE, 
                  cols = c("lightgrey", "red", "blue"), reduction = "umap") 

plot(p3[[3]])


### Check sexual stage marker expression

sextual.marker <- gsub('_', '-', Sextual.Stage$Bdiv)
p3 <- FeaturePlot(object = S.O.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = sextual.marker ,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p3)



###
p3 <- DotPlot(all.samples.integrated, features = my.gene, 
              cols = c("blue", "red", "green", "yellow"),
              dot.scale = 8) + RotatedAxis()

plot(p3)

p4 <- VlnPlot(object = all.samples.integrated, features = my.gene)
plot(p4)



# Differential Expression
# Identify global markers independent of cell cycle phase. Are there any global Cold Shock regulators?
# Two class comparison: fc = ident.1/ident.2

DefaultAssay(all.samples.integrated) <- "RNA"
Idents(all.samples.integrated) <- "spp"
objs <- unique(all.samples.integrated@meta.data$spp)
ref.obj <- objs[1]
quesy.obj <- objs[2:length(objs)]

global.shock.markers <- FindMarkers(all.samples.integrated, ident.1 = quesy.obj,  ident.2  = ref.obj, verbose = T)
global.shock.markers$genes <- rownames(global.shock.markers)
global.shock.markers$GeneID <- gsub('-', '_', global.shock.markers$genes)
global.shock.markers.sig <- global.shock.markers %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5))
## Percentages are calclulated on Variable features, not all genes
global.stats <- global.shock.markers.sig %>% 
  mutate(up.reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(up.reg) %>% summarise(num.DEGs = n(), percent = n() / nrow(global.shock.markers)) 

print(global.stats)

## Top up-regulated markers
global.shock.markers.top <- global.shock.markers.sig %>% dplyr::filter(avg_log2FC > 0) %>% top_n(4, abs(avg_log2FC))
global.shock.markers.top


p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = global.shock.markers.top$genes,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



p3 <- DotPlot(all.samples.integrated, features = global.shock.markers.top$genes, 
              cols = c("blue", "red", "green", "yellow"),
              dot.scale = 8) + RotatedAxis()

plot(p3)

VlnPlot(object = all.samples.integrated, features = global.shock.markers.top$genes)

tmp <- left_join(global.shock.markers.sig, prod.desc, by = 'GeneID')
write.xlsx(tmp, '../Output/scClockOut/global_shock_markers.xlsx')


### Differential expression analysis
## Cell cycle phase specific


all.spp.list <- SplitObject(all.samples.integrated, split.by = "spp")

## Re-normalize splitted data
for (i in 1:length(all.spp.list)) {
  Idents(all.spp.list[[i]]) <- 'phase.cond'
  all.spp.list[[i]] <- NormalizeData(all.spp.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}


GM <- getCellCyclePhaseMarkers(all.spp.list)


## For each cluster, pool global markers.
pooled_markers <- bind_rows(GM$all.markers.list.sig) %>% group_by(glob.clust) %>%
  summarise(glob.markers = list(unique(gene)), num.markers = length(unique(gene)))

objs <- unique(all.samples.integrated@meta.data$spp)
ref.obj <- objs[1]
quesy.obj <- objs[2:length(objs)]
clusters <- unique(all.samples.integrated@meta.data$seurat_clusters)
ref.clusts <- data.frame(ref = paste(ref.obj, clusters, sep = '_'))
ref.clusts$cluster <- gsub('.*_', '', ref.clusts$ref)
ref.clusts$dummy <- 1
query.clusts <- data.frame(query = paste(rep(quesy.obj, each = length(clusters)), clusters, sep = '_'))
query.clusts$dummy <- 1
query.clusts$cluster <- gsub('.*_', '', query.clusts$query)

contrasts <- full_join(ref.clusts, query.clusts, by = 'dummy') %>% 
  dplyr::filter(cluster.x == cluster.y) %>%
  transmute(ref = ref, query = query)

contrasts$cluster <- gsub('.*_', '', contrasts$ref)

#ident.1 case, ident.2 is control
Idents(all.samples.integrated) <- 'phase.cond'
matched.DEGs <- mclapply(split(contrasts, seq(nrow(contrasts))), function(x){
  tmp <- FindMarkers(all.samples.integrated, ident.1 = x$query, ident.2 = x$ref, verbose = T)
  ind <- rownames(tmp) %in% unlist(pooled_markers$glob.markers[which(pooled_markers$glob.clust == x$cluster)]) 
  tmp$ref <- x$ref
  tmp$query <- x$query
  tmp$cluster <- x$cluster
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp[ind, ])
  #return(tmp)
})


matched.DEGs <- bind_rows(matched.DEGs)

matched.DEGs.sig <- matched.DEGs %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5)) %>% arrange(desc(abs(avg_log2FC)))

matched.DEGs.stats <- matched.DEGs.sig %>%
  mutate(up.reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(query, up.reg) %>% summarise(num.DEGs = n()) 

print(matched.DEGs.stats)

matched.DEGs.stats$cluster <- as.factor(as.numeric(gsub('.*_', '', matched.DEGs.stats$query)))


matched.DEGs.top <- matched.DEGs.sig %>% dplyr::filter(avg_log2FC > 0) %>% group_by(cluster) %>% top_n(1, abs(avg_log2FC))



p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = matched.DEGs.top$gene,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



VlnPlot(object = all.samples.integrated, features = matched.DEGs.top$gene)


tmp <- left_join(matched.DEGs.sig, prod.desc, by = 'GeneID')
write.xlsx(tmp, '../Output/scClockOut/matched_shock_markers.xlsx')



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
prod.desc    <- read.xlsx('../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49   <- read.xlsx('../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx')
GT1.Pberghei <- read.xlsx('../Input/compScBdTgPb/Orthologs/rec_GT1.vs.P.bergei.xlsx')
GT1.bdiv     <- read.xlsx("../Input/compScBdTgPb/Orthologs/rec_GT1.vs.B.divergence.xlsx")

GT1.PB <- GT1.Pberghei %>% dplyr::select('query_id', contains("subject_id"))
colnames(GT1.PB) <- c('TGGT1', 'PBer')
GT1.BD <- GT1.bdiv %>% dplyr::select('query_id', contains('subject_id'))
colnames(GT1.BD) <- c('TGGT1', 'BDiv')

ME49.GT1.PB.BD <- inner_join(inner_join(TGGT1_ME49, GT1.PB, by = 'TGGT1'), GT1.BD, by = 'TGGT1')


## New Toxo Data
toxo.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.intra.expr.csv"

## PBer data
pb.file.csv  <- "../Input/MalariaCellAtlas/Expression_Matrices/10X/pb10xIDC/pb10xIDC_counts.csv"
pb.pheno.csv <- "../Input/MalariaCellAtlas/Expression_Matrices/10X/pb10xIDC/pb10xIDC_pheno.csv"


## BDiv_human
bdiv.file.csv  <- "../Input/scRNAseqBdiv/bdiv.expr.csv"


## Get expression matrix of orthologous genes
getExpr <- function(in.file, ME49.GT1.PB.BD, col.ind){
  file.counts <- read.csv(in.file)
  genes <- file.counts$X
  ind <- which(genes %in% ME49.GT1.PB.BD[,col.ind])
  file.counts <- file.counts[ind, ]
  genes <- genes[ind]
  genes <- ME49.GT1.PB.BD$TGGT1[match(genes, ME49.GT1.PB.BD[,col.ind])]
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  return(expr)
}

colnames(ME49.GT1.PB.BD)
toxo.counts <- getExpr(toxo.file.csv, ME49.GT1.PB.BD, 2)
pb.counts <- getExpr(pb.file.csv, ME49.GT1.PB.BD, 3)
bd.counts <- getExpr(bdiv.file.csv, ME49.GT1.PB.BD, 4)





## Generate Seurate objects
## individual Seurat objects
S.O.toxo <- CreateSeuratObject(counts = toxo.counts)
S.O.toxo$orig.ident <- 'Toxo'
S.O.toxo$spp <- 'Toxo'

S.O.pb <- CreateSeuratObject(counts = pb.counts)
S.O.pb$orig.ident <- 'Pberg'
S.O.pb$spp <- 'Pberg'

S.O.bd <- CreateSeuratObject(counts = bd.counts)
S.O.bd$orig.ident <- 'BDiv'
S.O.bd$spp <- 'BDiv'


S.O.list <- list(Toxo = S.O.toxo, PBerg = S.O.pb, BDiv = S.O.bd)

S.O.list <- mclapply(S.O.list, function(S.O){
  set.seed(100)
  S.O <- subset(x = S.O, downsample = 4000)
  
}, mc.cores = num.cores)

## Get common orthologous genes present in all data sets
comm.genes <- lapply(S.O.list, function(S.O){
  genes <- rownames(S.O@assays$RNA@data)
})

comm.genes <- Reduce(intersect, comm.genes)

S.O.list <- lapply(S.O.list, function(S.O){
  S.O <- subset(S.O, features = comm.genes)
})



## Quality Control

## Merge the data for qulity control
alldata <- merge(S.O.list[[1]], S.O.list[2:3], add.cell.ids=c("Toxo","PBerg","BDiv"))


feats <- c("nFeature_RNA","nCount_RNA")
VlnPlot(alldata, group.by= "orig.ident", features = feats, pt.size = 0.1,ncol = 4) + NoLegend()
FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

selected_c <- WhichCells(alldata, expression = nFeature_RNA > 100)
selected_f <- rownames(alldata)[ Matrix::rowSums(alldata) > 5]

## Filter the data
alldata <- subset(alldata, features=selected_f, cells=selected_c)
dim(alldata)


high.det <- WhichCells(alldata, expression = nFeature_RNA > 1200)

alldata <- subset(alldata, cells=setdiff(WhichCells(alldata),high.det))

ncol(alldata)



feats <- c("nFeature_RNA","nCount_RNA")
cowplot::plot_grid(ncol = 1,
                   VlnPlot(alldata, group.by= "orig.ident", features = feats, pt.size = 0.1,ncol = 4) + NoLegend())

Idents(alldata) <- 'spp'
p <- FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p


## Getting phenotypes

## For Toxo, transfer the labels from Boothroyed data
S.O.tg <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.RData')

## split the data, process each, transfer the lables
S.Os <- SplitObject(alldata, split.by = 'spp')
S.Os[[1]] <- prep_S.O(S.Os[[1]])
anchors <- FindTransferAnchors(reference = S.O.tg, query = S.Os[[1]], dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = S.O.tg@meta.data$phase,dims = 1:30)
predictions$phase <- predictions$predicted.id
S.Os[[1]] <- AddMetaData(object = S.Os[[1]], metadata = predictions)


## Reading Plasmodium Phenotypes
pb.pheno <- read.csv(pb.pheno.csv)
pb.pheno  <- pb.pheno  %>% 
  mutate(phase = case_when(absclust == 0 ~ "TrpE",
                           absclust == 1 ~ "TrpM",
                           absclust == 2 ~ "RngL",
                           absclust == 3 ~ "TrpL",
                           absclust == 4 ~ "SchE",
                           absclust == 5 ~ "SchL",
                           absclust == 6 ~ "RngE",
                           absclust == 7 ~ "SchM"))


S.Os[[2]] <- prep_S.O(S.Os[[2]])
rownames(pb.pheno) <- paste('PBerg_', pb.pheno$X, sep = '')
S.Os[[2]] <- AddMetaData(object = S.Os[[2]], metadata = pb.pheno)


## BDiv_human, no pheno available. Use Seurat clusters
S.Os[[3]] <- prep_S.O(S.Os[[3]])

S.Os[[3]]@meta.data$phase <- S.Os[[3]]@meta.data$seurat_clusters

## Remerge the data with new labels (metadata)
alldata.lab <- merge(S.Os[[1]], S.Os[2:3])


### Anchoring data to intra
S.O.list <- SplitObject(alldata.lab, split.by = 'spp')
S.O.list <- lapply(S.O.list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = S.O.list)
spps <- names(S.Os)
reference_dataset <- 1
anchors2 <- FindIntegrationAnchors(object.list = S.O.list, 
                                  anchor.features = features, 
                                  dims = 1:30)
                                  #reference = reference_dataset)

S.O.integrated <- IntegrateData(anchorset = anchors2, dims = 1:30)
# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"
S.O.integrated <- ScaleData(object = S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, features = VariableFeatures(object = S.O.integrated))
S.O.integrated <- FindNeighbors(S.O.integrated, dims = 1:10, reduction = 'pca')
S.O.integrated <- FindClusters(S.O.integrated, resolution = res)
S.O.integrated <- RunUMAP(S.O.integrated, dims = 1:13)



pcaMataData.integrated <- getPcaMetaData(S.O.integrated)


p1  <- ggplot(pcaMataData.integrated, aes(x= PC_1,y=-PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  #scale_color_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'crk2' = 'darkslateblue', 'ark3' = 'darkolivegreen4')) +
  #scale_fill_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'crk2' = 'darkslateblue', 'ark3' = 'darkolivegreen4')) +
  #scale_color_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +
  #scale_fill_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('UMAP2') + xlab('UMAP1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p1)

          
           

Idents(S.O.integrated) <- 'phase'
p <- DimPlot(S.O.integrated, reduction = "pca", 
             pt.size = 1,
             #shape.by='spp',
             split.by = 'spp',
             label = TRUE, label.size = 6) + #NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)


ggsave(filename="../Output/compScBdTgPb/figs/Tg_Pberg_Bdiv_pca_integrated_split.pdf", 
       plot=p,
       width = 12, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)




###### Bdiv Toxo infered phases.
Idents(S.O.integrated) <- 'spp'
S.O.toxo.bdiv <- subset(S.O.integrated, ident = c('Toxo', 'BDiv'))
Idents(S.O.toxo.bdiv) <- 'phase'
p <- DimPlot(S.O.toxo.bdiv, reduction = "pca", 
             pt.size = 1,
             #shape.by='spp',
             split.by = 'spp',
             label = TRUE, label.size = 6) + #NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)

Idents(S.O.toxo.bdiv) <- 'spp'
S.O.list.toxo.bdiv <- SplitObject(S.O.toxo.bdiv, split.by = 'spp')
anchors <- FindTransferAnchors(reference = S.O.list.toxo.bdiv[[1]], query = S.O.list.toxo.bdiv[[2]], dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = S.O.list.toxo.bdiv[[1]]@meta.data$phase,dims = 1:30)
predictions$phase <- predictions$predicted.id
S.O.list.toxo.bdiv[[2]] <- AddMetaData(object = S.O.list.toxo.bdiv[[2]], metadata = predictions)

Idents(S.O.list.toxo.bdiv[[2]]) <- 'phase'
p <- DimPlot(S.O.list.toxo.bdiv[[2]], reduction = "pca", 
             pt.size = 1,
             #shape.by='spp',
             split.by = 'spp',
             label = TRUE, label.size = 6) + #NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)


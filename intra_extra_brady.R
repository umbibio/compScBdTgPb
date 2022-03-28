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

## Count files
intra.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.intra.expr.csv"
extra.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.extra.expr.csv"
#crk2.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.crk2.expr.csv"
#ark3.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.ark3.expr.csv"

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
extra.counts <- getExpr(extra.file.csv, TGGT1_ME49)



abs.path <- "~/work/ToxoPlasmaGondiiR/Input/SingleCell/"
ID.Orthologs <- read.xlsx('~/work/ToxoPlasmaGondiiR/Input/ID_convert/convertIDs.xlsx')
SW3.meta.rds <- readRDS(file = paste(abs.path, "SW3_meta_data.rds",  sep = ""))
SW3.sparse.rds <- readRDS(file = paste(abs.path, "SW3_sparse_expression.rds",  sep = ""))

M <- as.data.frame(as.matrix(SW3.sparse.rds))
Meta <- as.data.frame(as.matrix(SW3.meta.rds))
Meta <- Meta %>% mutate(Sample = rownames(Meta))

## Convert IDs to GT1
M.GT1 <- M[!is.na(match(rownames(M), ID.Orthologs$TGME49ID)),]
rownames(M.GT1) <- ID.Orthologs$TGGT1ID[match(rownames(M.GT1), ID.Orthologs$TGME49ID)]

Meta.72 <- Meta %>% dplyr::filter(Timepoint == 72)
M.72 <- M.GT1 %>% dplyr::select(colnames(M)[colnames(M.GT1) %in% Meta.72$Sample])

#Meta.24 <- Meta %>% dplyr::filter(Timepoint == 24)

#M.72.KO.D3 <- M.72 %>% dplyr::select(contains('KO.D3'))
#M.72.KO.pH <- M.72 %>% dplyr::select(contains('KO'),-contains('D3'))
#M.72.WT.D3 <- M.72 %>% dplyr::select(contains('WT.D3'))
M.72.WT.pH <- M.72 %>% dplyr::select(contains('WT.pH')) ## Bradyzoites
M.72.WT.D3 <- M.72 %>% dplyr::select(contains('WT.D3')) ## Tachyzoites
#M.72.WT    <- M.72 %>% dplyr::select(contains('WT')) ## contains both stressed and unstressed
#M.72.D3    <- M.72 %>% dplyr::select(contains('D3'))


S.O.intra <- CreateSeuratObject(counts = intra.counts, min.cells = 3, min.features = 100)
S.O.intra$orig.ident <- 'intra'
VlnPlot(S.O.intra, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.intra, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.intra  <- subset(S.O.intra , subset = nFeature_RNA > 100 & nFeature_RNA < 1200 )

S.O.extra <- CreateSeuratObject(counts = extra.counts, min.cells = 3, min.features = 100)
S.O.extra$orig.ident <- 'extra'
VlnPlot(S.O.extra, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.extra, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.extra  <- subset(S.O.extra , subset = nFeature_RNA > 100 & nFeature_RNA < 1300 )


S.O.pru.WT <- CreateSeuratObject(counts = M.72.WT.D3, min.cells = 3, min.features = 100)
S.O.pru.WT$orig.ident <- 'pru.WT' 
VlnPlot(S.O.pru.WT, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.pru.WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.pru.WT  <- subset(S.O.pru.WT, subset = nFeature_RNA > 100 & nFeature_RNA < 1500 )

S.O.pru.WT.pH <- CreateSeuratObject(counts = M.72.WT.pH, min.cells = 3, min.features = 100)
S.O.pru.WT.pH$orig.ident <- 'pru.pH'
VlnPlot(S.O.pru.WT.pH, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.pru.WT.pH, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.pru.WT.pH  <- subset(S.O.pru.WT.pH, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 )


#### Quick test on Brady data
Brady.Pru <- prep_S.O(S.O.pru.WT.pH)
Brady.Pru@meta.data$spp <- 'Brady'

Tachy.Pru <- prep_S.O(S.O.pru.WT)
Tachy.Pru@meta.data$spp <- 'Tachy'

DimPlot(Tachy.Pru, reduction = 'pca')
DimPlot(Brady.Pru, reduction = 'umap')

SAG1.ID <- 'TGGT1-233460'
BAG1.ID <- "TGGT1-259020"


p <- FeaturePlot(Tachy.Pru, features = c(SAG1.ID, BAG1.ID))
plot(p)

p <- FeaturePlot(Brady.Pru, features = c(SAG1.ID, BAG1.ID))
plot(p)


Tachy.Brady <- merge(Tachy.Pru, Brady.Pru, add.cell.ids=c("Tachy","Brady"))
Tachy.Brady.list <- SplitObject(Tachy.Brady, split.by = "spp")
Tachy.Brady.features <- SelectIntegrationFeatures(object.list = Tachy.Brady.list)
Tachy.Brady.anchors <- FindIntegrationAnchors(object.list = Tachy.Brady.list, anchor.features = Tachy.Brady.features)
Tachy.Brady.integrated <- IntegrateData(anchorset = Tachy.Brady.anchors)
DefaultAssay(Tachy.Brady.integrated) <- "integrated"
Tachy.Brady.integrated <- ScaleData(Tachy.Brady.integrated, verbose = FALSE)
Tachy.Brady.integrated <- RunPCA(Tachy.Brady.integrated, npcs = 30, verbose = FALSE)
Tachy.Brady.integrated <- RunUMAP(Tachy.Brady.integrated, reduction = "pca", dims = 1:30)
Tachy.Brady.integrated <- FindNeighbors(Tachy.Brady.integrated, reduction = "pca", dims = 1:30)
Tachy.Brady.integrated <- FindClusters(Tachy.Brady.integrated, resolution = 0.2)
Idents(Tachy.Brady.integrated) <- 'spp'
DimPlot(Tachy.Brady.integrated, reduction = 'umap', label = T)

DefaultAssay(Tachy.Brady.integrated) <- "RNA"
p <- FeaturePlot(Tachy.Brady.integrated, features = c(SAG1.ID, BAG1.ID))
plot(p)
#######
#S.O.list <- list(intra = S.O.intra, extra = S.O.extra, pru.WT = S.O.pru.WT, pru.pH = S.O.pru.WT.pH)
#S.O.list <- list(intra = S.O.intra, extra = S.O.extra, pru.WT = S.O.pru.WT)
S.O.list <- list(intra = S.O.intra, extra = S.O.extra, pru.pH = S.O.pru.WT.pH)

S.O.list <- mclapply(S.O.list, function(S.O){
  set.seed(100)
  S.O <- subset(x = S.O, downsample = 6000)
  
}, mc.cores = num.cores)

## Get common orthologous genes present in all data sets
comm.genes <- lapply(S.O.list, function(S.O){
  genes <- rownames(S.O@assays$RNA@data)
})

comm.genes <- Reduce(intersect, comm.genes)

S.O.list <- lapply(S.O.list, function(S.O){
  S.O <- subset(S.O, features = comm.genes)
})




## Merge the data for qulity control
#alldata <- merge(S.O.list[[1]], S.O.list[2:4], add.cell.ids=c("intra","extra","pru.WT","pru.pH"))
#alldata <- merge(S.O.list[[1]], S.O.list[2:3], add.cell.ids=c("intra","extra","pru.WT"))
alldata <- merge(S.O.list[[1]], S.O.list[2:3], add.cell.ids=c("intra","extra","pru.pH"))

alldata$spp <- alldata$orig.ident



## Boothroyed data
S.O.tg <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.RData')

## split the data, process each, transfer the lables
S.Os <- SplitObject(alldata, split.by = 'spp')
S.Os <- mclapply(S.Os, function(S.O){
  S.O <- prep_S.O(S.O)
  anchors <- FindTransferAnchors(reference = S.O.tg, query = S.O, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.tg@meta.data$phase,dims = 1:30)
  predictions$phase <- predictions$predicted.id
  #predictions$phase[which(predictions$prediction.score.max < 0.7)] <- 'NA'
  S.O <- AddMetaData(object = S.O, metadata = predictions)
  return(S.O)
}, mc.cores = num.cores)

spps <- names(S.Os)

S.Os <- lapply(1:length(S.Os), function(i){
  S.Os[[i]]@meta.data$spp <- spps[i]
  S.Os[[i]]
})


#saveRDS(S.Os, '../Input/compScBdTgPb/RData/S.O.intra_extra_pruWT_pruPH_lables_not_anchored.RData')


## Test plots
Idents(S.Os[[3]]) <- 'phase'
p <- DimPlot(S.Os[[3]], reduction = "umap", 
             #group.by = "cell", 
             #split.by = 'spp',
             pt.size = 1,
             #shape.by='spp',
             label = F, label.size = 4) + #NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)


SAG1.ID <- 'TGGT1-233460'
BAG1.ID <- "TGGT1-259020"
BFD1.ID <- "TGGT1-200385"
new.feature <- "TGGT1-311100"

p <- FeaturePlot(S.Os[[3]], features = c(SAG1.ID, BAG1.ID, BFD1.ID, new.feature))
plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/brady_sag1_bag1_expr.pdf", 
       plot=p,
       width = 8, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Extract bradyzoites
#BAG1_expression = GetAssayData(object = S.Os[[3]], assay = "RNA", slot = "data")["TGGT1-259020",]
#pos_ids = names(which(BAG1_expression>0))
#neg_ids = names(which(CD14_expression==0))

## Merge the data for qulity control
#alldata.phase <- merge(S.Os[[1]], S.Os[2:4], add.cell.ids=c("intra","extra","pru.WT","pru.pH"))
alldata.phase <- merge(S.Os[[1]], S.Os[2:3], add.cell.ids=c("intra","extra","pru.pH"))

## Integrate the data
S.O.list <- SplitObject(alldata.phase, split.by = "spp")
S.O.list <- lapply(X = S.O.list, FUN = function(x) {
  ## Extract the count data
  
  ## extract the count data from each as.matrix(S.O.list[[1]][["RNA"]]@data)
  ## Replace genes with Bdiv orthologous when needed
  ## recreate the new Seurat object.
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = S.O.list)
anchors <- FindIntegrationAnchors(object.list = S.O.list, anchor.features = features)
S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
S.O.integrated <- ScaleData(S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, npcs = 30, verbose = FALSE)
S.O.integrated <- RunUMAP(S.O.integrated, reduction = "pca", dims = 1:30)
S.O.integrated <- FindNeighbors(S.O.integrated, reduction = "pca", dims = 1:30)
S.O.integrated <- FindClusters(S.O.integrated, resolution = 0.2)

S.O.integrated$phase.spp <- paste(S.O.integrated@meta.data$spp, S.O.integrated@meta.data$phase, sep = "_")
Idents(S.O.integrated) <- "phase.spp"

DimPlot(S.O.integrated, reduction = 'umap', label = T, split.by = 'spp') + NoLegend()

## Non-anchord version
S.O.no.anchored <- merge(S.O.list[[1]], S.O.list[2:3], add.cell.ids=c("intra","extra","pru.pH"), merge.data = T)
S.O.no.anchored <- FindVariableFeatures(S.O.no.anchored, selection.method = "vst", nfeatures = 3000)
S.O.no.anchored <- ScaleData(S.O.no.anchored)
S.O.no.anchored <- RunPCA(S.O.no.anchored, features = VariableFeatures(object = S.O.no.anchored))
S.O.no.anchored <- FindNeighbors(S.O.no.anchored, dims = 1:30, reduction = 'pca')
S.O.no.anchored <- FindClusters(S.O.no.anchored, resolution = 0.2)
S.O.no.anchored <- RunUMAP(S.O.no.anchored, dims = 1:30)

S.O.no.anchored$phase.spp <- paste(S.O.no.anchored@meta.data$spp, S.O.no.anchored@meta.data$phase, sep = "_")
Idents(S.O.no.anchored) <- "spp"

DimPlot(S.O.no.anchored, reduction = 'umap', group.by = 'phase', label = T, label.size = 6) + 
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
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = 'None',
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold"))


## Extracting the data for control over PCA directions
getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp, phase = S.O@meta.data$phase)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}



pcaMataData.alldata <- getPcaMetaData(S.O.no.anchored)

pcaMataData.alldata$phase <- factor(pcaMataData.alldata$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))

p1  <- ggplot(pcaMataData.alldata, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase, 
  ), #color = 'blue', 
  alpha = 0.9,
  shape=21, size = 0.3)+ 
  #scale_color_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  #scale_fill_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  
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
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = c(0.1, 0.2),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))


plot(p1)


ggsave(filename="../Output/compScBdTgPb/figs/intra_extra_brady_umap.pdf", 
       plot=p1,
       width = 8, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

## Filter S/M/C out of Extra and Brady
S.O.no.anchored.filter <-  S.O.no.anchored
S.O.no.anchored.filter@meta.data$keep <- 'yes'
ind1 <- S.O.no.anchored.filter@meta.data$spp == 'pru.pH' & !(S.O.no.anchored.filter@meta.data$phase %in% c("G1.a", "G1.b"))
ind2 <- S.O.no.anchored.filter@meta.data$spp == 'extra' & !(S.O.no.anchored.filter@meta.data$phase %in% c("G1.a", "G1.b"))

S.O.no.anchored.filter@meta.data$keep[ind1] <- 'no'
S.O.no.anchored.filter@meta.data$keep[ind2] <- 'no'
Idents(S.O.no.anchored.filter) <- 'keep'
S.O.no.anchored.filter <- subset(S.O.no.anchored.filter, idents = 'yes')


pcaMataData.alldata.filter <- getPcaMetaData(S.O.no.anchored.filter)

pcaMataData.alldata.filter$phase <- factor(pcaMataData.alldata.filter$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))

p1  <- ggplot(pcaMataData.alldata.filter, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase, 
  ), #color = 'blue', 
  alpha = 0.9,
  shape=21, size = 0.3)+ 
  #scale_color_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  #scale_fill_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  
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
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = c(0.1, 0.2),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))


plot(p1)




ggsave(filename="../Output/compScBdTgPb/figs/intra_extra_brady_umap_filter.pdf", 
       plot=p1,
       width = 8, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

## Look at filtered integrated data for marker analysis

S.O.integrated.filter <-  S.O.integrated
S.O.integrated.filter@meta.data$keep <- 'yes'
ind1 <- S.O.integrated.filter@meta.data$spp == 'pru.pH' & !(S.O.integrated.filter@meta.data$phase %in% c("G1.a", "G1.b"))
ind2 <- S.O.integrated.filter@meta.data$spp == 'extra' & !(S.O.integrated.filter@meta.data$phase %in% c("G1.a", "G1.b"))

S.O.integrated.filter@meta.data$keep[ind1] <- 'no'
S.O.integrated.filter@meta.data$keep[ind2] <- 'no'
Idents(S.O.integrated.filter) <- 'keep'
S.O.integrated.filter <- subset(S.O.integrated.filter, idents = 'yes')

Idents(S.O.integrated.filter) <- 'phase'
p <- DimPlot(S.O.integrated.filter, reduction = "umap", 
             #group.by = "cell", 
             split.by = 'spp',
             pt.size = 1,
             #shape.by='spp',
             label = T, label.size = 4) + 
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
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = 'None',
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold"))



plot(p)

pcaMataData.S.O.integrated.filter <- getPcaMetaData(S.O.integrated.filter)

pcaMataData.S.O.integrated.filter$phase <- factor(pcaMataData.S.O.integrated.filter$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))

p1  <- ggplot(pcaMataData.S.O.integrated.filter, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = spp,
    color = spp, 
  ), #color = 'blue', 
  alpha = 0.9,
  shape=21, size = 0.3)+ 
  #scale_color_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  #scale_fill_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  
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
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = c(0.1, 0.2),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))


plot(p1)


ggsave(filename="../Output/compScBdTgPb/figs/intra_extra_brady_umap_filter.pdf", 
       plot=p1,
       width = 8, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


DefaultAssay(S.O.integrated.filter) <- 'RNA'
Idents(S.O.integrated.filter) <- 'spp'

markers <- FindAllMarkers(object = S.O.integrated.filter, only.pos = T, min.pct = 0.3) 
colnames(markers)[colnames(markers) == 'cluster'] <- 'spp'
markers$GeneID <- gsub('-', '_', markers$gene)
markers.top <- markers %>% group_by(spp) %>% slice_max(n = 1, order_by = avg_log2FC)

markers.sig <- markers %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.05) 
Idents(S.O.integrated.filter) <- 'phase'
p <- FeaturePlot(object = S.O.integrated.filter, features = markers.top$gene, label = T,repel = T,
                 #shape.by  = 'spp',
                 split.by = 'spp',
                 cols = c("grey", "blue"), reduction = "umap")
plot(p)

markers.sig <- left_join(markers.sig, prod.desc, by = 'GeneID') 
markers.sig <- markers.sig %>% arrange(spp, desc(avg_log2FC))

write.xlsx(markers.sig , '../Output/compScBdTgPb/tables/extra_vs_intra_vs_brady_global_markers.xlsx')



markers.stat <- markers.sig %>% group_by(spp) %>% summarise(num.deg = n())

# bar plot 
p <-  ggplot(markers.stat, aes(x=spp, y=num.deg)) +
  geom_bar(stat="identity", fill = "steelblue")+
  geom_text(aes(label=num.deg), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Intra", "Extra", "Brady")) +
  theme(axis.text.x = element_text(face="bold", size=14, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=14, angle=0, color = 'black'),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/intra_extra_brady_total_degs.pdf", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Intra RH vs Intra Pru
#Intra.vs.Pru.markers <- FindMarkers(S.O.integrated.G1, ident.1 = "intra", ident.2 = "pru.WT", min.pct = 0.3)
#Intra.vs.Pru.markers$gene <- rownames(Intra.vs.Pru.markers)
#Intra.vs.Pru.markers$GeneID <- gsub('-', '_', Intra.vs.Pru.markers$gene)
#Intra.vs.Pru.markers.sig <- Intra.vs.Pru.markers %>% dplyr::filter(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05) 
#Intra.vs.Pru.markers.sig <- left_join(Intra.vs.Pru.markers.sig, prod.desc, by = 'GeneID') 
#write.xlsx(Intra.vs.Pru.markers.sig , '../Output/compScBdTgPb/tables/Intra_vs_PruWT_markers.xlsx')

## Intra RH vs Brady
Idents(S.O.integrated.filter) <- 'spp'
Intra.vs.Brady.markers <- FindMarkers(S.O.integrated.filter, ident.1 = "intra", ident.2 = "pru.pH", min.pct = 0.3)
Intra.vs.Brady.markers$gene <- rownames(Intra.vs.Brady.markers)
Intra.vs.Brady.markers$GeneID <- gsub('-', '_', Intra.vs.Brady.markers$gene)
Intra.vs.Brady.markers.sig <- Intra.vs.Brady.markers %>% dplyr::filter(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05) 
Intra.vs.Brady.markers.sig <- left_join(Intra.vs.Brady.markers.sig, prod.desc, by = 'GeneID') 
write.xlsx(Intra.vs.Brady.markers.sig , '../Output/compScBdTgPb/tables/Intra_vs_brady_markers.xlsx')


## Intra RH vs Extra RH
Intra.vs.Extra.markers <- FindMarkers(S.O.integrated.filter, ident.1 = "intra", ident.2 = "extra", min.pct = 0.3)
Intra.vs.Extra.markers$gene <- rownames(Intra.vs.Extra.markers)
Intra.vs.Extra.markers$GeneID <- gsub('-', '_', Intra.vs.Extra.markers$gene)
Intra.vs.Extra.markers.sig <- Intra.vs.Extra.markers %>% dplyr::filter(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05) 
Intra.vs.Extra.markers.sig <- left_join(Intra.vs.Extra.markers.sig, prod.desc, by = 'GeneID') 
write.xlsx(Intra.vs.Brady.markers.sig , '../Output/compScBdTgPb/tables/Intra_vs_Extra_markers.xlsx')

## Brady vs Extra
#Brady.vs.Extra.markers <- FindMarkers(S.O.integrated.G1, ident.1 = "pru.pH", ident.2 = "extra", min.pct = 0.3)
#Brady.vs.Extra.markers$gene <- rownames(Brady.vs.Extra.markers)
#Brady.vs.Extra.markers$GeneID <- gsub('-', '_', Brady.vs.Extra.markers$gene)
#Brady.vs.Extra.markers.sig <- Brady.vs.Extra.markers %>% dplyr::filter(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05) 
#Brady.vs.Extra.markerssig <- left_join(Brady.vs.Extra.markers.sig, prod.desc, by = 'GeneID') 
#write.xlsx(Brady.vs.Extra.markers.sig , '../Output/compScBdTgPb/tables/Brady_vs_Extra_markers.xlsx')


## Identify RH vs Brady specific gene sets (remove RH vs Pru WT)
#tmp <- inner_join(Intra.vs.Pru.markers.sig, Intra.vs.Brady.markers.sig, by = 'GeneID')
#tmp$dir <- sign(tmp$avg_log2FC.x * tmp$avg_log2FC.y)
#RH.Pru.genes <- tmp$GeneID[tmp$dir == 1]

#Intra.vs.Brady.markers.sig.specific <- Intra.vs.Brady.markers.sig %>% dplyr::filter(!(GeneID %in% RH.Pru.genes))
#write.xlsx(Intra.vs.Brady.markers.sig.specific , '../Output/compScBdTgPb/tables/Intra_vs_brady_markers_specific.xlsx')


Intra.vs.Brady.Extra <- full_join(Intra.vs.Extra.markers.sig, Intra.vs.Brady.markers.sig, by = 'GeneID')
Intra.vs.Brady.Extra$dir <- sign(Intra.vs.Brady.Extra$avg_log2FC.x * Intra.vs.Brady.Extra$avg_log2FC.y)
Intra.vs.Brady.Extra <- Intra.vs.Brady.Extra %>% arrange(desc(dir), avg_log2FC.x)
colnames(Intra.vs.Brady.Extra) <- gsub('y', 'Brady', gsub('x', 'Extra', colnames(Intra.vs.Brady.Extra)))
Intra.vs.Brady.Extra$avg_log2FC.Extra <- -Intra.vs.Brady.Extra$avg_log2FC.Extra
Intra.vs.Brady.Extra$avg_log2FC.Brady <- -Intra.vs.Brady.Extra$avg_log2FC.Brady
write.xlsx(Intra.vs.Brady.Extra , '../Output/compScBdTgPb/tables/Intra_vs_brady_extra_markers.xlsx')


## Intra RH vs Extra RH (NO FILTERS)
Idents(S.O.integrated) <- 'spp'
Intra.vs.Extra.no.phase.filt.markers <- FindMarkers(S.O.integrated, ident.1 = "intra", ident.2 = "extra", min.pct = 0.3)
Intra.vs.Extra.no.phase.filt.markers$gene <- rownames(Intra.vs.Extra.no.phase.filt.markers)
Intra.vs.Extra.no.phase.filt.markers$GeneID <- gsub('-', '_', Intra.vs.Extra.no.phase.filt.markers$gene)
Intra.vs.Extra.no.phase.filt.markers.sig <- Intra.vs.Extra.no.phase.filt.markers %>% dplyr::filter(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05) 
Intra.vs.Extra.no.phase.filt.markers.sig <- left_join(Intra.vs.Extra.no.phase.filt.markers.sig, prod.desc, by = 'GeneID') 

#### Feature plots
my.genes <- c(
  'TGGT1-208370',
  'TGGT1-214980',
  'TGGT1-295935',
  'TGGT1-212300',
  'TGGT1-311100',
  'TGGT1-226380',
  'TGGT1-316250',
  'TGGT1-273320',
  'TGGT1-221590',
  'TGGT1-248990',
  'TGGT1-288360',
  'TGGT1-225540',
  'TGGT1-313910',
  'TGGT1-320490',
  'TGGT1-204340',
  'TGGT1-321530',
  'TGGT1-311230',
  'TGGT1-289540',
  'TGGT1-247030',
  'TGGT1-233380',
  'TGGT1-294870',
  'TGGT1-215910',
  'TGGT1-248740',
  'TGGT1-236650',
  'TGGT1-289600',
  'TGGT1-244560',
  'TGGT1-273760',
  'TGGT1-204050')

gene.names <- c(
  'MHC',
  'hypothetical protein',
  'KRUF',
  'hypothetical protein',
  'zinc finger (CCCH type)',
  'hypothetical protein',
  'hypothetical protein',
  'hypothetical protein',
  'dual specificity phosphatase',
  'hypothetical protein',
  'tryptophanyl-tRNA synthetase (TrpRS2)',
  'hypothetical protein',
  'RNA recognition motif 2 protein',
  'N-acyl-phosphatidylethanolamine-hydrolyzing phospholipase D family protein',
  'hypothetical protein',
  'cathepsin CPL',
  'hypothetical protein',
  'hypothetical protein',
  'hypothetical protein',
  'hypothetical protein',
  'universal stress family protein',
  'hypothetical protein',
  'hypothetical protein',
  'DEAD (Asp-Glu-Ala-Asp) box polypeptide 17',
  'HSP29',
  'HSP90',
  'HSP70',
  'SUB1')
S.O.no.anchored.filter$spp2 <- S.O.no.anchored.filter$spp
S.O.no.anchored.filter$spp2[S.O.no.anchored.filter$spp == 'intra'] <- 'Intra'
S.O.no.anchored.filter$spp2[S.O.no.anchored.filter$spp == 'extra'] <- 'Extra'
S.O.no.anchored.filter$spp2[S.O.no.anchored.filter$spp == 'pru.pH'] <- 'Brady'
S.O.no.anchored.filter$spp2 <- factor(S.O.no.anchored.filter$spp2, levels = c('Intra', 'Extra', 'Brady'))
Idents(S.O.no.anchored.filter) <- 'spp2'
DefaultAssay(S.O.no.anchored.filter) <- 'RNA'

i <- 1
i <- 5
i <- 9
i <- 13
i <- 21

FeaturePlot(object = S.O.no.anchored.filter, features = 'TGGT1-312100', label = T,repel = T)
p <- FeaturePlot(object = S.O.no.anchored.filter, features = my.genes[i], label = F,repel = T,label.size = 6, 
                 #shape.by  = 'spp',
                 #split.by = 'spp',
                 cols = c("grey", "blue"), reduction = "umap") + 
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
  #facet_wrap(spp~.) + 
  ggtitle(paste(gene.names[i], gsub('-', '_', my.genes[i]),sep = ':')) +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = c(0.1, 0.25),
        #legend.position = 'None',
        legend.title = element_text(colour="black", size=6, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=6, 
                                   face="bold"))



plot(p)

p <- grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave(filename="../Output/compScBdTgPb/figs/MHC_expr.pdf", 
       plot=p,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)





##### Feature Plots for important markers on unfiltered data
S.O.no.anchored$spp2 <- S.O.no.anchored$spp
S.O.no.anchored$spp2[S.O.no.anchored$spp == 'intra'] <- 'Intra'
S.O.no.anchored$spp2[S.O.no.anchored$spp == 'extra'] <- 'Extra'
S.O.no.anchored$spp2[S.O.no.anchored$spp == 'pru.pH'] <- 'Brady'
S.O.no.anchored$spp2 <- factor(S.O.no.anchored$spp2, levels = c('Intra', 'Extra', 'Brady'))
Idents(S.O.no.anchored) <- 'spp2'
DefaultAssay(S.O.no.anchored) <- 'RNA'



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
i <- 1
p <- FeaturePlot(object = S.O.no.anchored, features = my.genes[i], label = F,repel = T,label.size = 6, 
                 #shape.by  = 'spp',
                 #split.by = 'spp',
                 cols = c("grey", "blue"), reduction = "umap") + 
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
  #facet_wrap(spp~.) + 
  ggtitle(paste(gene.names[i], gsub('-', '_', my.genes[i]),sep = ':')) +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = c(0.1, 0.25),
        #legend.position = 'None',
        legend.title = element_text(colour="black", size=6, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=6, 
                                   face="bold"))



plot(p)

### Violin of important markers

show.lab2 <- c('TGGT1_259020', # BAG1
               'TGGT1_306620', # AP2IX-9
               'TGGT1_200385', # BFD1
               'TGGT1_233460' # SAG1
)

show.lab2.names <- c('BAG1', 'AP2IX-9', 'BFD1', 'SAG1')
VlnPlot(subset(S.O.no.anchored.filter, idents = c('Brady', 'Intra')), features = gsub('_', '-', show.lab))
p <- VlnPlot(S.O.no.anchored.filter, features = gsub('_', '-', show.lab2))


p <- lapply(1:4, function(i){
  p[[i]] + 
    ylab('Expr') + xlab('') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    #facet_wrap(spp~.) + 
    ggtitle(paste(show.lab2.names[i], gsub('-', '_', show.lab2[i]),sep = ':')) +
    theme(
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) + 
    theme(#legend.position = c(0.1, 0.25),
      legend.position = 'None',
      legend.title = element_text(colour="black", size=6, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=6, 
                                 face="bold"))
})

p.all <- grid.arrange(p[[1]], p[[4]], p[[2]], p[[3]], ncol = 2)

ggsave(filename="../Output/compScBdTgPb/figs/AP2_BFD1_BAG1_SAG1_violin.pdf",
       plot=p.all,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## Quantify the overlap
library(ComplexHeatmap)
library(UpSetR)
library(ggVennDiagram)

## Using Intra as base
brady.up <- Intra.vs.Brady.markers.sig %>% dplyr::filter(avg_log2FC < 0) %>% arrange(avg_log2FC)
brady.down <- Intra.vs.Brady.markers.sig %>% dplyr::filter(avg_log2FC > 0) %>% arrange(avg_log2FC)
Extra.up <- Intra.vs.Extra.markers.sig %>% dplyr::filter(avg_log2FC < 0) %>% arrange(avg_log2FC)
Extra.down <- Intra.vs.Extra.markers.sig %>% dplyr::filter(avg_log2FC > 0) %>% arrange(avg_log2FC)

Brady.Extra <- list(Brady.up = brady.up$GeneID, Brady.down = brady.down$GeneID, 
                    Extra.up = Extra.up$GeneID, Extra.down = Extra.down$GeneID)

g <- ggVennDiagram(Brady.Extra, label_size = 5, set_size = 6,label = "count", label_alpha = 0,label_color = 'white',
                   category.names = c("Brady (Up)", "Brady (Down)", "Extra (Up)", "Extra (Down)"),
                   set_color = c("brady.up" = "firebrick","brady.down" ="darkorchid3", 'Extra.up' = 'darkslateblue', 
                                 'Extra.down' = 'darkolivegreen4')) +  
  scale_color_manual(values = c("firebrick","darkorchid3", 'darkslateblue', 'darkolivegreen4')) + 
  scale_x_continuous(expand = expansion(mult = .15)) + theme(legend.position = "none")

plot(g)


ggsave(filename="../Output/compScBdTgPb/figs/Brady_Extra_shared_degs.pdf", 
       plot=g,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

#######

## Enrichment Analysis

in.dir <- '../Input/GO/brady_tachy_extra/'
all.files <- list.files(in.dir)

all.clust.items <- list()
for(f in all.files){
  nn <- gsub('\\.tsv', '', f)
  spp <- strsplit(nn, split = '_')[[1]][1]
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp$spp <- spp
  all.clust.items <- c(all.clust.items, list(tmp))
}


all.clust.items <- do.call(rbind, all.clust.items)

filtered.Go <- all.clust.items %>% arrange(spp) %>% distinct() %>%
  group_by(spp) %>% mutate(rank = row_number()) %>%
  dplyr::filter(Benjamini < 0.1 & rank < 25) %>% 
  arrange(spp) 

filtered.Go$Name <- gsub('process', 'proc.', filtered.Go$Name)
filtered.Go$spp <- factor(filtered.Go$spp, levels = c('Intra', 'Extra', 'Brady'))
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level=unique(filtered.Go$Name))
## Category of contrasts
p <- ggplot(filtered.Go, aes(x = spp, y = Name)) + 
  geom_point(aes(colour = "red", size = -log(Benjamini))) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) 

plot(p)


ggsave(filename="../Output/compScBdTgPb/figs/intra_extra_brady_GO.pdf", 
       plot=p,
       width = 6, height = 10, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


### Intra vs Brady

Intra.vs.Brady.markers.all <- FindMarkers(S.O.integrated.filter, ident.1 = "intra", ident.2 = "pru.pH", min.pct = 0)
Intra.vs.Brady.markers.all$gene <- rownames(Intra.vs.Brady.markers.all)
Intra.vs.Brady.markers.all$GeneID <- gsub('-', '_', Intra.vs.Brady.markers.all$gene)
volc.df <- left_join(Intra.vs.Brady.markers.all, prod.desc, by = 'GeneID')
volc.df$sig <- ifelse(abs(volc.df$avg_log2FC) > 0.58 & volc.df$p_val_adj < 0.01, 'yes', 'no')
volc.df$geneLable <- volc.df$ProductDescription
volc.df$geneLable <- gsub('microneme protein ', '', gsub('rhoptry protein ', '', 
                          gsub('SRS29B', 'SAG1', 
                               gsub('AP2 domain transcription factor ', '',
                                    gsub('SAG-related sequence ', '' , 
                                         gsub('bradyzoite antigen ', '', volc.df$geneLable))))))
volc.df$geneLable[which(volc.df$GeneID == 'TGGT1_200385')] <- 'BFD1'
show.lab <- c('TGGT1_259020', # BAG1
              'TGGT1_280570', # SRS35A
              'TGGT1_264660', # SRS44
              'TGGT1_202020', # DnAK-TPR
              'TGGT1_306620', # AP2IX-9
              'TGGT1_269010', #AP2VIII-7
              'TGGT1_240900', #AP2VI-2
              'TGGT1_280460', # AP2VIIa-2
              'TGGT1_200385', # BFD1
              'TGGT1_233460', # SAG1
              'TGGT1_271050', # SRS34A
              'TGGT1_309590', #ROP1
              'TGGT1_291890' # MIC1
              )

# volc.df <- volc.df %>% group_by(contrast) %>% 
#   mutate(GeneLable = ifelse(-log(qval) > 20  & abs(log_fc) > 1.7 , "True", "False"))
volc.df <- volc.df %>% mutate(GeneLable = ifelse(GeneID %in%  show.lab, "True", "False"))


p <- ggplot(volc.df) + geom_point(aes(avg_log2FC, -log10(p_val_adj + 1e-301),color=sig))+ 
  geom_text_repel(aes(avg_log2FC, -log10(p_val_adj + 1e-301)), size = 5, fontface = "bold",
                  label = ifelse(volc.df$GeneLable == "True", as.character(volc.df$geneLable),""), 
                  box.padding = unit(0.45, "lines"),
                  hjust=1,
                  segment.angle = 120,
                  nudge_x = -1.5, 
                  nudge_y = 1.5,
                  segment.size = 0.4) + 
  theme(legend.title=element_blank(),text = element_text(size=20))+ 
  scale_color_manual(values = c("yes" = "red", "no" = "black"))+
  theme_bw()+
  theme(
    #plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    #legend.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 18, angle = 90, vjust = 0.45, color = "black"),  
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, color = "black", face = "bold"),
    strip.text = element_text(color = "black", size = 18, face = "bold"),
    strip.background = element_rect(fill = "white"))+
  theme(legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  xlab("log2 FC")

p


ggsave(filename="../Output/coldShock/figs/volcano.pdf",
       plot=p,
       width = 10, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

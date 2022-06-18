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

set.seed(100)
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



S.O.intra <- CreateSeuratObject(counts = intra.counts, min.cells = 3, min.features = 100)
S.O.intra$orig.ident <- 'intra'
S.O.intra$spp <- 'intra'
#VlnPlot(S.O.intra, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O.intra, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.intra  <- subset(S.O.intra , subset = nFeature_RNA > 100 & nFeature_RNA < 1200 )
S.O.intra  <- subset(S.O.intra, downsample = 1500)

S.O.extra <- CreateSeuratObject(counts = extra.counts, min.cells = 3, min.features = 100)
S.O.extra$orig.ident <- 'extra'
S.O.extra$spp <- 'extra'
#VlnPlot(S.O.extra, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O.extra, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.extra  <- subset(S.O.extra , subset = nFeature_RNA > 100 & nFeature_RNA < 1300 )
S.O.extra  <- subset(S.O.extra, downsample = 1500)

S.O.intra.extra <- merge(S.O.intra, S.O.extra, add.cell.ids=c("intra","extra"))

## Boothroyed RH data
S.O.tg <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.RData')

## split the data, process each, transfer the labels
S.Os <- SplitObject(S.O.intra.extra, split.by = 'spp')
S.Os <- mclapply(S.Os, function(S.O){
  S.O <- prep_S.O(S.O)
  anchors <- FindTransferAnchors(reference = S.O.tg, query = S.O, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.tg@meta.data$phase,dims = 1:30)
  predictions$phase <- predictions$predicted.id
  #predictions$phase[which(predictions$prediction.score.max < 0.7)] <- 'NA'
  S.O <- AddMetaData(object = S.O, metadata = predictions)
  return(S.O)
}, mc.cores = num.cores)

DimPlot(S.Os[[1]], reduction = 'umap')


saveRDS(S.Os, '../Input/compScBdTgPb/RData/S.O.intra.extra.labs_list.rds')

## Boothroyd Pru data
## Clusters 7 & 8 correspond to cluster P6
S.O.pru.list <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.pru.boothroyd.list2.rds')


### Merge all data
S.O.list <- list(intra = S.Os[[1]], extra = S.Os[[2]], 
                 pru.d0 = S.O.pru.list$d0, pru.d3 = S.O.pru.list$d3,
                 pru.d5 = S.O.pru.list$d5, pru.d7 = S.O.pru.list$d7)


S.O.list <- mclapply(S.O.list, function(S.O){
  set.seed(100)
  S.O <- subset(x = S.O, downsample = 4000)
  
}, mc.cores = num.cores)


## Get common orthologous genes present in all data sets
# comm.genes <- lapply(S.O.list, function(S.O){
#   genes <- rownames(S.O@assays$RNA@data)
# })
# 
# comm.genes <- Reduce(intersect, comm.genes)
# 
# S.O.list <- lapply(S.O.list, function(S.O){
#   S.O <- subset(S.O, features = comm.genes)
# })




## Merge the data for qulity control
alldata <- merge(S.O.list[[1]], S.O.list[2:6], add.cell.ids=c("intra","extra","pru.d0","pru.d3",
                                                              "pru.d5", 'pru.d7'))

alldata$spp <- alldata$orig.ident
alldata$match <- ifelse(alldata$spp %in% c('intra', 'extra'), alldata$spp, alldata$cluster)

saveRDS(alldata, '../Input/compScBdTgPb/RData/S.O.intra_extra_pru_d0_d7_lables_list_boothroyd.RData')

## Integrate the data
S.O.list <- SplitObject(alldata, split.by = "spp")
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

### Adding better spp names
S.O.integrated@meta.data$spp2 <- ifelse(S.O.integrated@meta.data$spp == 'intra', 'RH.intra',
                                             ifelse(S.O.integrated@meta.data$spp == 'extra', 'RH.extra',
                                                    ifelse(S.O.integrated@meta.data$spp == 'pru.d0', 'Pru.intra',
                                                           'Pru.brady')))
S.O.integrated@meta.data$spp2 <- factor(S.O.integrated@meta.data$spp2, 
                                        levels = c('RH.intra', 'RH.extra', 'Pru.intra', 'Pru.brady'))
S.O.integrated$phase.spp2 <- paste(S.O.integrated@meta.data$spp2, S.O.integrated@meta.data$phase, sep = "_")
Idents(S.O.integrated) <- "phase.spp2"

Idents(S.O.integrated) <- "spp2"
DimPlot(S.O.integrated, reduction = 'umap', label = T,label.size = 6) + NoLegend()

## matched by looking at ids from here: http://st-atlas.org/
## and matching the metadata renamed_clusters
S.O.integrated$spp3 <- ifelse(S.O.integrated$match == 1, 'P3', 
                              ifelse(S.O.integrated$match == 2, 'P4',
                                     ifelse(S.O.integrated$match == 3, 'P4', 
                                            ifelse(S.O.integrated$match == 4, 'P5', 
                                                   ifelse(S.O.integrated$match == 5, 'P5',
                                                          ifelse(S.O.integrated$match == 6, 'P5',
                                                                 ifelse(S.O.integrated$match == 7, 'P6',
                                                                        ifelse(S.O.integrated$match == 8, 'P6',
                                                                               ifelse(S.O.integrated$match == 9, 'P1',
                                                                                      ifelse(S.O.integrated$match == 'intra',
                                                                                             'RH.intra', 'RH.extra'))))))))))



Idents(S.O.integrated) <- "spp3"
DimPlot(S.O.integrated, reduction = 'umap', label = T) + NoLegend()

saveRDS(S.O.integrated, '../Input/compScBdTgPb/RData/S.O.intra.extra.pru_d0_d7.integrated_boothroyd.rds')
## Extracting the data for control over PCA directions
## Extracting the data for control over PCA directions
getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp2, phase = S.O@meta.data$phase, 
                          match = S.O@meta.data$spp3)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}

pcaMataData.alldata <- getPcaMetaData(S.O.integrated)
pcaMataData.alldata$phase <- factor(pcaMataData.alldata$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
pcaMataData.alldata$spp <- factor(pcaMataData.alldata$spp, levels = c('RH.intra', 'RH.extra', 'Pru.intra', 'Pru.brady'))
pcaMataData.alldata$match <- factor(pcaMataData.alldata$match, levels = c('RH.intra', 'RH.extra', 
                                                              'P1', 'P2', 'P3', 'P4', 'P5', 'P6'))



p2  <- ggplot(pcaMataData.alldata, aes(x= UMAP_1,y=UMAP_2)) +
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
  facet_wrap(spp~.) + 
  theme(
    strip.text.x = element_text(
      size = 14,  face = "bold.italic"
    ),
    strip.text.y = element_text(
      size = 14, face = "bold.italic"
    )
  ) +
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = "None",
    #legend.position = c(0.14, 0.14),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))


plot(p2)


ggsave(filename="../Output/compScBdTgPb/figs/boothroyd/intra_extra_brady_umap_boothroyd_split.pdf", 
       plot=p2,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Feature plots
DefaultAssay(S.O.integrated) <- 'RNA'

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
              'TGGT1-233460',
              'TGGT1-228170')


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
                "SAG1",
                'GRA44')


#i <- 3
i <- 14
p <- FeaturePlot(object = S.O.integrated, features = my.genes[i], label = F,repel = T,label.size = 6,
                 #shape.by  = 'spp',
                 split.by = 'spp2',
                 ncol = 1,
                 cols = c("grey", "blue"), reduction = "umap")

p.dat <- lapply(1:4, function(i) p[[i]]$data)
p.dat <- bind_rows(p.dat)
colnames(p.dat) <- c('UMAP1', 'UMAP2', 'spp', 'expr')
g <- ggplot(data = p.dat, aes(x = UMAP1, y = UMAP2, color = expr, fill = expr)) + 
  geom_point(aes(#fill = lable.prob,
    #fill = expr,
    #color = 'gray66', 
  ), #color = 'blue', 
  alpha = 0.9,
  shape=21, size = 0.4)+ 
  theme_bw(base_size = 14) +
  scale_color_gradient2(midpoint = 0.01, low = "white", mid = "gray",
                        high = "blue", space = "Lab" ) + 
  scale_fill_gradient2(midpoint = 0.01, low = "white", mid = "gray",
                        high = "blue", space = "Lab" ) + 
  ylab('UMAP2') + xlab('UMAP1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  facet_wrap(spp~.) + 
  theme(
    strip.text.x = element_text(
      size = 14,  face = "bold.italic"
    ),
    strip.text.y = element_text(
      size = 14, face = "bold.italic"
    )
  ) +
  
  ggtitle(paste(gene.names[i], gsub('-', '_', my.genes[i]),sep = ':')) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = "None",
        #legend.position = c(0.88, 0.17),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=2)))


 
plot(g)


ggsave(filename='../Output/compScBdTgPb/figs/boothroyd/AP2IX-10_expr.pdf',
       plot=g,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)






#####

prep_un_anchored <- function(S.O){
  S.O <- FindVariableFeatures(S.O, selection.method = "vst", nfeatures = 3000)
  S.O <- ScaleData(S.O)
  S.O <- RunPCA(S.O, features = VariableFeatures(object = S.O))
  S.O <- FindNeighbors(S.O, dims = 1:30, reduction = 'pca')
  S.O <- FindClusters(S.O, resolution = 0.2)
  S.O <- RunUMAP(S.O, dims = 1:30)
  S.O
}
## Non-anchord version
S.O.no.anchored <- prep_un_anchored(alldata)
S.O.no.anchored$phase.spp <- paste(S.O.no.anchored@meta.data$spp, S.O.no.anchored@meta.data$phase, sep = "_")
Idents(S.O.no.anchored) <- "spp"

DimPlot(S.O.no.anchored, reduction = 'umap')
## Filter out S/M/C from Extra
S.O.no.anchored2 <- merge(S.O.list[[1]], S.O.list[2:4], add.cell.ids=c("intra","extra","pru.WT", "pru.pH"), merge.data = T)
S.O.no.anchored2@meta.data$keep <- 'yes'
ind <- S.O.no.anchored2@meta.data$spp == 'extra' & !(S.O.no.anchored2@meta.data$phase %in% c("G1.a", "G1.b"))
S.O.no.anchored2@meta.data$keep[ind] <- 'no'
Idents(S.O.no.anchored2) <- 'keep'
S.O.no.anchored2 <- subset(S.O.no.anchored2, idents = 'yes')
S.O.no.anchored2 <- prep_un_anchored(S.O.no.anchored2)
S.O.no.anchored2$phase.spp <- paste(S.O.no.anchored2@meta.data$spp, S.O.no.anchored2@meta.data$phase, sep = "_")
Idents(S.O.no.anchored2) <- "spp"


DimPlot(S.O.no.anchored2, reduction = 'umap',  label = T, label.size = 6) + 
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


pcaMataData.alldata <- getPcaMetaData(S.O.no.anchored2)

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
  theme(legend.position = c(0.1, 0.8),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))


plot(p1)


ggsave(filename="../Output/compScBdTgPb/figs/intra_extra_brady_umap2.pdf", 
       plot=p1,
       width = 8, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

## Filter S/M/C out of All for DEG analysis
S.O.integrated.filter <-  S.O.integrated
Idents(S.O.integrated.filter) <- 'phase'
S.O.integrated.filter <- subset(S.O.integrated.filter, idents = c('G1.a', 'G1.b'))
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

pcaMataData.S.O.integrated.filter$phase <- factor(pcaMataData.S.O.integrated.filter$phase, levels = c('G1.a', 'G1.b'))

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


## DEG analysis
DefaultAssay(S.O.integrated.filter) <- 'RNA'

## Combine intra and pruWT for DEG analysis vs. Brady
S.O.integrated.filter@meta.data$stage = ifelse(S.O.integrated.filter@meta.data$spp %in% c('intra', 'pru.WT'), 'Tachy',
                                               ifelse(S.O.integrated.filter@meta.data$spp == 'pru.pH', 'Brady', 'Extra'))
Idents(S.O.integrated.filter) <- 'stage'
Global.Markers <- FindAllMarkers(object = S.O.integrated.filter, only.pos = T, min.pct = 0.3) 
colnames(Global.Markers)[colnames(Global.Markers) == 'cluster'] <- 'stage'
Global.Markers$GeneID <- gsub('-', '_', Global.Markers$gene)
Global.Markers.top <- Global.Markers %>% group_by(stage) %>% slice_max(n = 1, order_by = avg_log2FC)

Global.Markers.sig <- Global.Markers %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.05) 
p <- FeaturePlot(object = S.O.integrated.filter, features = Global.Markers.top$gene, label = T,repel = T,
                 #shape.by  = 'spp',
                 split.by = 'stage',
                 cols = c("grey", "blue"), reduction = "umap")
plot(p)

Global.Markers.sig <- left_join(Global.Markers.sig, prod.desc, by = 'GeneID') 
Global.Markers.sig <- Global.Markers.sig %>% arrange(stage, desc(avg_log2FC))

write.xlsx(Global.Markers.sig, '../Output/compScBdTgPb/tables/Tachy_Brady_Extra_global_markers.xlsx')



Global.Markers.stat <- Global.Markers.sig %>% group_by(stage) %>% summarise(num.deg = n())

# bar plot 
p <-  ggplot(Global.Markers.stat, aes(x=stage, y=num.deg)) +
  geom_bar(stat="identity", fill = "steelblue")+
  geom_text(aes(label=num.deg), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Tachy", "Extra", "Brady")) +
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

ggsave(filename="../Output/compScBdTgPb/figs/Tachy_Extra_Brady_total_degs.pdf", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Pairwise comparisons
## Tachy vs Brady
Idents(S.O.integrated.filter) <- 'stage'
Brady.vs.Tachy.markers <- FindMarkers(S.O.integrated.filter, ident.1 = "Brady", ident.2 = "Tachy", min.pct = 0.3)
Brady.vs.Tachy.markers$gene <- rownames(Brady.vs.Tachy.markers)
Brady.vs.Tachy.markers$GeneID <- gsub('-', '_', Brady.vs.Tachy.markers$gene)
Brady.vs.Tachy.markers.sig <- Brady.vs.Tachy.markers %>% dplyr::filter(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05) 
Brady.vs.Tachy.markers.sig <- left_join(Brady.vs.Tachy.markers.sig, prod.desc, by = 'GeneID') 
write.xlsx(Brady.vs.Tachy.markers.sig , '../Output/compScBdTgPb/tables/Brady_vs_Tachy_markers.xlsx')


## Tachy RH vs Extra RH 
Idents(S.O.integrated.filter) <- 'spp'
Extra.vs.Intra.RH.markers <- FindMarkers(S.O.integrated.filter, ident.1 = "extra", ident.2 = "intra", min.pct = 0.3)
Extra.vs.Intra.RH.markers$gene <- rownames(Extra.vs.Intra.RH.markers)
Extra.vs.Intra.RH.markers$GeneID <- gsub('-', '_', Extra.vs.Intra.RH.markers$gene)
Extra.vs.Intra.RH.markers.sig <- Extra.vs.Intra.RH.markers %>% dplyr::filter(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05) 
Extra.vs.Intra.RH.markers.sig <- left_join(Extra.vs.Intra.RH.markers.sig, prod.desc, by = 'GeneID') 
write.xlsx(Extra.vs.Intra.RH.markers.sig , '../Output/compScBdTgPb/tables/Extra_vs_Intra_RH_markers.xlsx')


## Shared Brady/Extra markers compared to Tachy
shared.Brady.Extra <- inner_join(Brady.vs.Tachy.markers.sig, Extra.vs.Intra.RH.markers.sig, by = 'GeneID')
shared.Brady.Extra$dir <- sign(shared.Brady.Extra$avg_log2FC.x * shared.Brady.Extra$avg_log2FC.y)
shared.Brady.Extra <- shared.Brady.Extra %>% arrange(desc(dir), avg_log2FC.x)
colnames(shared.Brady.Extra) <- gsub('x', 'Brady', gsub('y', 'Extra', colnames(shared.Brady.Extra)))
write.xlsx(shared.Brady.Extra, '../Output/compScBdTgPb/tables/shared_brady_extra_markers.xlsx')



## Venn Comparison
library(ComplexHeatmap)
library(UpSetR)
library(ggVennDiagram)

## Using Intra as base
brady.up <- Brady.vs.Tachy.markers.sig %>% dplyr::filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
brady.down <- Brady.vs.Tachy.markers.sig %>% dplyr::filter(avg_log2FC < 0) %>% arrange(avg_log2FC)
Extra.up <- Extra.vs.Intra.RH.markers.sig %>% dplyr::filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
Extra.down <- Extra.vs.Intra.RH.markers.sig %>% dplyr::filter(avg_log2FC < 0) %>% arrange(avg_log2FC)

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

i <- 1
i <- 5
i <- 9
i <- 13
i <- 21

p <- FeaturePlot(object = S.O.no.anchored2, features = my.genes[i], label = F,repel = T,label.size = 6, 
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
  theme(legend.position = c(0.1, 0.75),
        #legend.position = 'None',
        legend.title = element_text(colour="black", size=6, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=6, 
                                   face="bold"))



plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/MHC_expr.pdf", 
       plot=p,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


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

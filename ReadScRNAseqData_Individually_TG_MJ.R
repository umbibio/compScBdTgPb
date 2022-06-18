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

## KO genes
Crk2 <- 'TGGT1-218220'
Ark3 <- 'TGGT1-203010'

## Count files
intra.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.intra.expr.csv"
extra.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.extra.expr.csv"
crk2.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.crk2.expr.csv"
ark3.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.ark3.expr.csv"

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
crk2.counts <- getExpr(crk2.file.csv, TGGT1_ME49)
ark3.counts <- getExpr(ark3.file.csv, TGGT1_ME49)


## individual Seurat objects
S.O.intra <- CreateSeuratObject(counts = intra.counts)
S.O.intra$orig.ident <- 'intra'
S.O.extra <- CreateSeuratObject(counts = extra.counts)
S.O.extra$orig.ident <- 'extra'
S.O.crk2  <- CreateSeuratObject(counts = crk2.counts)
S.O.crk2$orig.ident <- 'crk2'
S.O.ark3  <- CreateSeuratObject(counts = ark3.counts)
S.O.ark3$orig.ident <- 'ark3'

S.O.list <- list(intra = S.O.intra, extra = S.O.extra, crk2 = S.O.crk2, ark3 = S.O.ark3)

S.O.list <- mclapply(S.O.list, function(S.O){
  set.seed(100)
  S.O <- subset(x = S.O, downsample = 6000)
  
}, mc.cores = num.cores)


## Merge the data for qulity control
alldata <- merge(S.O.list[[1]], S.O.list[2:4], add.cell.ids=c("intra","extra","crk2","ark3"))
alldata$spp <- alldata$orig.ident

# as.data.frame(alldata@assays$RNA@counts[1:10,1:2])
# head(alldata@meta.data,10)

feats <- c("nFeature_RNA","nCount_RNA")
VlnPlot(alldata, group.by= "orig.ident", features = feats, pt.size = 0.1,ncol = 4) + NoLegend()
FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

selected_c <- WhichCells(alldata, expression = nFeature_RNA > 100)
selected_f <- rownames(alldata)[ Matrix::rowSums(alldata) > 5]
selected_f <- unique(c(selected_f, Crk2, Ark3)) ## Make sure Crk2 and Ark 3 are selected
                
## Filter the data
alldata <- subset(alldata, features=selected_f, cells=selected_c)
dim(alldata)


high.det <- WhichCells(alldata, expression = nFeature_RNA > 1200)

alldata <- subset(alldata, cells=setdiff(WhichCells(alldata),high.det))

ncol(alldata)


# rel_expression <- t( t(as.matrix(alldata@assays$RNA@counts)) / Matrix::colSums(as.matrix(alldata@assays$RNA@counts))) * 100
# most_expressed <- sort(Matrix::rowSums( rel_expression ),T)[20:1] / ncol(alldata)
# 
# par(mfrow=c(1,2),mar=c(4,6,1,1))
# boxplot( as.matrix(t(rel_expression[names(most_expressed),])),cex=.1, las=1, xlab="% total count per cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)


feats <- c("nFeature_RNA","nCount_RNA")
cowplot::plot_grid(ncol = 1,
                   VlnPlot(alldata, group.by= "orig.ident", features = feats, pt.size = 0.1,ncol = 4) + NoLegend())

Idents(alldata) <- 'spp'
p <- FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p


## Individually process the data and transfer labels

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

saveRDS(S.Os, '../Input/compScBdTgPb/RData/S.O.intra_extra_crk2_ark3_lables_not_anchored.RData')
saveRDS(S.Os[[1]], '../Input/compScBdTgPb/RData/S.O.toxo_MJ_lables.RData')

## For updating the APP
saveRDS(S.Os[[1]], '../scExpressionProfiler/S_O_toxo_MJ_labels.rds')

## Test plots
Idents(S.Os[[1]]) <- 'spp'
p <- DimPlot(S.Os[[1]], reduction = "umap", 
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

## Extracting the data for control over PCA directions
getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), phase = S.O@meta.data$phase, 
                          predicted.id = S.O@meta.data$predicted.id, 
                          lable.prob = S.O@meta.data$prediction.score.max,
                          spp = S.O@meta.data$spp)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}


## Plot lable transfered individual S.Os
pcaMetaData <- lapply(S.Os, getPcaMetaData)
pcaMetaData <- bind_rows(pcaMetaData)
pcaMetaData$predicted.id <- factor(pcaMetaData$predicted.id, 
                                   levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
pcaMetaData$phase <- factor(pcaMetaData$phase, 
                                   levels = c('G1.a', 'G1.b', 'S', 'M', 'C', 'NA'))

pcaMetaData$spp <- factor(pcaMetaData$spp, 
                            levels = c('intra', 'extra', 'crk2', 'ark3'))







##### Phase proportions

## Proportions
tmp <- pcaMetaData

stats <- tmp %>% group_by(spp) %>% mutate(total.cells = n()) %>% 
  ungroup() %>% group_by(spp, phase) %>% summarise(counts = n(), perc = n()/total.cells[1]) 

aes(x=spp, y=perc, color = phase, group = phase)
p3 <- ggplot(data=stats, aes(x=spp, y=perc, fill = phase)) +
  geom_bar(stat="identity", position=position_stack(), width=0.8) +
  geom_text(aes(label=round(perc, 2)), vjust=1.2,  color="black", 
            size=5, fontface="bold", position=position_stack(), angle=0)+
  #scale_fill_manual(values = c("G1" = "#ff6c67","S/M" ='#a2a700', 'C' = '#00c377')) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=16, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=16, angle=0)) +
  theme(
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) #+  coord_flip()

plot(p3)



ggsave(filename="../Output/compScBdTgPb/figs/phase_proportions_intra_extra_crk2_ark3.pdf",
       plot=p3,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


###### Flipping plots

getPlot <- function(pcaMetaData, reduction = 'umap', flip1 = 1, flip2 = 1){
  if(reduction == 'umap'){
    p  <- ggplot(pcaMetaData, aes(x=flip1 * UMAP_1,y=flip2 * UMAP_2)) +
      geom_point(aes(#fill = lable.prob,
        fill = predicted.id,
        color = predicted.id
      ), #color = 'blue', 
      shape=21, size = 1)+ 
      scale_color_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +
      scale_fill_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +
                                
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
  }else{
    p  <- ggplot(pcaMetaData, aes(x=flip1 * PC_1,y=flip2 * PC_2)) +
      geom_point(aes(#fill = lable.prob,
        fill = predicted.id,
        color = predicted.id
      ), #color = 'blue', 
      shape=21, size = 1)+ 
      scale_color_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +
      scale_fill_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +
      
      theme_bw(base_size = 14) +
      theme(legend.position = "right") +
      #scale_fill_gradientn(colours = viridis::inferno(10)) +
      #scale_fill_gradientn(colours = col_range(10)) +
      #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
      #scale_fill_brewer(palette = "BuPu") +
      ylab('PC2') + xlab('PC1') +
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
  }
  
  return(p)
}



pcaMetaData.intra <- pcaMetaData %>% dplyr::filter(spp == 'intra')
pcaMetaData.extra <- pcaMetaData %>% dplyr::filter(spp == 'extra')
pcaMetaData.extra$predicted.id <- factor(pcaMetaData.extra$predicted.id, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
pcaMetaData.crk2  <- pcaMetaData %>% dplyr::filter(spp == 'crk2')
pcaMetaData.ark3  <- pcaMetaData %>% dplyr::filter(spp == 'ark3')


p1 <- getPlot(pcaMetaData.intra,  reduction = 'umap', flip1 = 1, flip2 = 1)
p2 <- getPlot(pcaMetaData.extra,  reduction = 'umap', flip1 = 1, flip2 = -1)
p3 <- getPlot(pcaMetaData.crk2,  reduction = 'umap', flip1 = 1, flip2 = 1)
p4 <- getPlot(pcaMetaData.ark3,  reduction = 'umap', flip1 = 1, flip2 = -1)
p <- grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave(filename="../Output/compScBdTgPb/figs/individually_processed_umap_label_transfer.pdf",
       plot=p,
       width = 12, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


p1 <- getPlot(pcaMetaData.intra,  reduction = 'pca', flip1 = 1, flip2 = -1)
p2 <- getPlot(pcaMetaData.extra,  reduction = 'pca', flip1 = 1, flip2 = -1)
p3 <- getPlot(pcaMetaData.crk2,  reduction = 'pca', flip1 = 1, flip2 = 1)
p4 <- getPlot(pcaMetaData.ark3,  reduction = 'pca', flip1 = 1, flip2 = 1)
p <- grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave(filename="../Output/compScBdTgPb/figs/individually_processed_pca_label_transfer.pdf",
       plot=p,
       width = 12, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


### Merging data
alldata.lab <- merge(S.Os[[1]], S.Os[2:4], add.cell.ids=c("intra","extra","crk2","ark3"))
alldata.lab$spp <- alldata.lab$orig.ident

alldata.lab <- prep_S.O(alldata.lab)

pcaMataData.alldata.lab <- getPcaMetaData(alldata.lab)


p1  <- ggplot(pcaMataData.alldata.lab, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = spp,
    color = spp
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'crk2' = 'darkslateblue', 'ark3' = 'darkolivegreen4')) +
  scale_fill_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'crk2' = 'darkslateblue', 'ark3' = 'darkolivegreen4')) +
  
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
  )


p2  <- ggplot(pcaMataData.alldata.lab, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = predicted.id,
    color = predicted.id
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +
  scale_fill_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +

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
  )


p <- grid.arrange(p1, p2, ncol = 1)

ggsave(filename="../Output/compScBdTgPb/figs/merged_transfered_labels.pdf",
       plot=p,
       width = 8, height = 12,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


### Anchoring data to intra
## shared PCA/UMAP

S.O.list <- lapply(X = S.Os, FUN = function(x) {
  ## Extract the count data
  
  ## extract the count data from each as.matrix(S.O.list[[1]][["RNA"]]@data)
  ## Replace genes with Bdiv orthologous when needed
  ## recreate the new Seurat object.
  DefaultAssay(x) <- 'RNA'
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = S.O.list)
spps <- names(S.Os)
reference_dataset <- 1
anchors <- FindIntegrationAnchors(object.list = S.O.list, 
                                  anchor.features = features, reference = reference_dataset)
S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"
S.O.integrated <- ScaleData(object = S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, features = VariableFeatures(object = S.O.integrated))
S.O.integrated <- FindNeighbors(S.O.integrated, dims = 1:10, reduction = 'pca')
S.O.integrated <- FindClusters(S.O.integrated, resolution = res)
S.O.integrated <- RunUMAP(S.O.integrated, dims = 1:13)


pcaMataData.integrated <- getPcaMetaData(S.O.integrated)


p1  <- ggplot(pcaMataData.integrated, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = predicted.id,
    color = predicted.id
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  #scale_color_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'crk2' = 'darkslateblue', 'ark3' = 'darkolivegreen4')) +
  #scale_fill_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'crk2' = 'darkslateblue', 'ark3' = 'darkolivegreen4')) +
  scale_color_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +
  scale_fill_manual(values = c("G1.a" = "firebrick","G1.b" ="darkorange2", 'S' = 'gold3', 'M' = 'darkolivegreen4', 'C' = 'darkorchid2')) +
  
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

## Global spp specific markers
DefaultAssay(S.O.integrated) <- 'RNA'
Idents(S.O.integrated) <- 'spp'
S.O.extra.intra <- subset(S.O.integrated, ident = c('intra', 'extra'))

DefaultAssay(S.O.extra.intra) <- 'RNA'
Idents(S.O.extra.intra) <- 'spp'

markers.extra.intra <- FindAllMarkers(object = S.O.extra.intra, only.pos = T, min.pct = 0.1) 
colnames(markers.extra.intra)[colnames(markers.extra.intra) == 'cluster'] <- 'spp'
markers.extra.intra$GeneID <- gsub('-', '_', markers.extra.intra$gene)
markers.extra.intra.top <- markers.extra.intra %>% group_by(spp) %>% slice_max(n = 1, order_by = avg_log2FC)

markers.extra.intra.sig <- markers.extra.intra %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.05) 
p <- FeaturePlot(object = S.O.extra.intra, features = 'TGGT1-215895', label = T,repel = T,
                 #shape.by  = 'spp',
                 split.by = 'spp',
                 cols = c("grey", "blue"), reduction = "pca")
plot(p)

markers.extra.intra.sig <- left_join(markers.extra.intra.sig, prod.desc, by = 'GeneID') 
markers.extra.intra.sig <- markers.extra.intra.sig %>% arrange(spp, desc(avg_log2FC))

write.xlsx(markers.extra.intra.sig , '../Output/compScBdTgPb/tables/extra_vs_intra_global_markers.xlsx')





## Global spp specific markers ALL SPPS
DefaultAssay(S.O.integrated) <- 'RNA'
Idents(S.O.integrated) <- 'spp'


markers.all <- FindAllMarkers(object = S.O.integrated, only.pos = T, min.pct = 0.15) 
colnames(markers.all)[colnames(markers.all) == 'cluster'] <- 'spp'
markers.all$GeneID <- gsub('-', '_', markers.all$gene)
markers.all.top <- markers.all %>% group_by(spp) %>% slice_max(n = 1, order_by = avg_log2FC)

markers.all.sig <- markers.all %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.05) 
p <- FeaturePlot(object = S.O.integrated, features = markers.all.top$gene, label = T,repel = T,
                 #shape.by  = 'spp',
                 split.by = 'spp',
                 cols = c("grey", "blue"), reduction = "pca")
plot(p)

markers.all.sig <- left_join(markers.all.sig, prod.desc, by = 'GeneID') 
markers.all.sig <- markers.all.sig %>% arrange(spp, desc(avg_log2FC))

write.xlsx(markers.all.sig , '../Output/compScBdTgPb/tables/all_spp_global_markers.xlsx')

#### Here we need to generate a bar plot of spp specific global markers total numbers.
markers.all.sig.stats <- markers.all.sig %>% group_by(spp) %>% summarise(num.deg = n())
markers.all.sig.stats$spp <- factor(markers.all.sig.stats$spp, levels = c('intra', 'extra', 'crk2', 'ark3'))
p <- ggplot(data=markers.all.sig.stats, aes(x=spp, y=num.deg)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() + 
  geom_text(aes(label=num.deg), vjust=1.6, color="black", size=4.5, fontface = 'bold')+
  ylab('upregulated genes') + xlab('spicies') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) 
  
 


plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/global_markers_deg_numbers.pdf",
       plot=p,
       width = 8, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## Marker analysis for infered phases
spps <- names(S.Os)
all.markers.list <- mclapply(1:length(spps), function(i){
  S.O <- S.Os[[i]]
  DefaultAssay(S.O) <- "RNA"
  Idents(S.O) <- 'predicted.id'
  markers <- FindAllMarkers(object = S.O, only.pos = TRUE)
  markers$GeneID = gsub('-', '_', markers$gene)
  colnames(markers)[colnames(markers) == 'cluster'] <- 'predicted.id'
  markers$spp <- spps[i]
  return(markers)
}, mc.cores = num.cores) 

all.markers <- bind_rows(all.markers.list)

all.markers.sig <- all.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
all.markers.sig.stats <- all.markers.sig %>% group_by(spp, predicted.id) %>% summarise(num.deg = n())

all.markers.sig.stats$spp <- factor(all.markers.sig.stats$spp, 
                                    levels = c('intra', 'extra', 'crk2', 'ark3'))
all.markers.sig.stats$predicted.id <- factor(all.markers.sig.stats$predicted.id, 
                                   levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))


p <- ggplot(data=all.markers.sig.stats, aes(x=predicted.id, y=num.deg)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.deg), vjust=1.6, color="black", size=5)+
  theme_bw() + 
  ylab('number of markers') + xlab('predicted phase') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  
  facet_wrap(spp~.)

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/individually_processed_degs_of_transfered_labels.pdf",
       plot=p,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## Marker similarity analysis
all.markers.sig.list <- all.markers.sig %>% ungroup() %>% group_by(spp, predicted.id) %>% summarise(genes = list(GeneID))
all.markers.sig.list$dummy <- 1
all.markers.sig.list <- full_join(all.markers.sig.list, all.markers.sig.list, by = 'dummy')

spp.ranks <- data.frame(spp = c('intra', 'extra', 'crk2', 'ark3'), rank = c(1,2,3,4))
all.markers.sig.list <- left_join(all.markers.sig.list, spp.ranks, by = c("spp.x" = "spp"))
all.markers.sig.list <- left_join(all.markers.sig.list, spp.ranks, by = c("spp.y" = "spp"))

all.markers.sig.list <- all.markers.sig.list %>% dplyr::filter(rank.x <= rank.y)
all.markers.sig.list <- all.markers.sig.list %>% rowwise() %>% 
  mutate(overlap = length(intersect(unlist(genes.x), unlist(genes.y))))

hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab')
all.markers.sig.list$spp.x <- factor(all.markers.sig.list$spp.x, 
                                     levels = c('intra', 'extra', 'crk2', 'ark3'))
all.markers.sig.list$spp.y <- factor(all.markers.sig.list$spp.y, 
                                     levels = c('intra', 'extra', 'crk2', 'ark3'))
all.markers.sig.list$predicted.id.x <- factor(all.markers.sig.list$predicted.id.x, 
                                             levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
all.markers.sig.list$predicted.id.y <- factor(all.markers.sig.list$predicted.id.y, 
                                              levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))


p <- ggplot(all.markers.sig.list, aes(x = predicted.id.x, y = predicted.id.y, fill = overlap)) + 
  geom_tile(aes(fill = overlap)) + 
  geom_text(aes(label = overlap), size = 4, fontface='bold') + 
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) + theme_minimal() + 
  ylab("") + xlab("") + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 0, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1,face = 'bold', color = 'black')) +
  
  scale_fill_gradientn(colours = hm.palette(10)) + 
  facet_grid(spp.y~spp.x,switch = "both") + 
  theme(
    strip.text.x = element_text(
      size = 12, color = "red", face = "bold.italic"
    ),
    strip.text.y = element_text(
      size = 12, color = "red", face = "bold.italic"
    )
  )
  
plot(p)  



ggsave(filename="../Output/compScBdTgPb/figs/cross_compare_degs_transfered_labels.pdf.pdf", plot=p,
       width = 8, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Marker expression plots
my.genes <- c(Crk2, Ark3)

Idents(alldata.lab) <- 'predicted.id'
p <- VlnPlot(alldata.lab, my.genes, group.by = 'predicted.id')
plot(p)

alldata.lab$predicted.id <- factor(alldata.lab$predicted.id, level = c('G1.a', 'G1.b', 'S', 'M', 'C'))
alldata.lab$spp <- factor(alldata.lab$spp, level = c('intra', 'extra', 'crk2', 'ark3'))
p1 <- DotPlot(alldata.lab, features = my.genes, group.by = 'spp') + RotatedAxis()
p2 <- DotPlot(alldata.lab, features = my.genes, group.by = 'predicted.id') + RotatedAxis()
p1+p2

p <- grid.arrange(p1, p2, ncol = 2)
ggsave(filename="../Output/compScBdTgPb/figs/crk2_ark3_expr.pdf", plot=p,
       width = 10, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

top.markers <- all.markers.sig %>% group_by(predicted.id, spp) %>% slice_max(n = 1, order_by = avg_log2FC)

Idents(S.Os[[1]]) <- 'predicted.id'
Idents(S.Os[[2]]) <- 'predicted.id'
Idents(S.Os[[3]]) <- 'predicted.id'
Idents(S.Os[[4]]) <- 'predicted.id'

p1 <- FeaturePlot(S.Os[[1]], features = Crk2, split.by = 'spp', reduction = 'umap', label = T)
p2 <- FeaturePlot(S.Os[[2]], features = Crk2, split.by = 'spp', reduction = 'umap', label = T)
p3 <- FeaturePlot(S.Os[[3]], features = Crk2, split.by = 'spp', reduction = 'umap', label = T)
p4 <- FeaturePlot(S.Os[[4]], features = Crk2, split.by = 'spp', reduction = 'umap', label = T)

p <- grid.arrange(p1,p2,p3,p4, ncol = 2)

ggsave(filename="../Output/compScBdTgPb/figs/ark3_expr.pdf", plot=p,
       width = 8, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


p1 <- FeaturePlot(S.Os[[1]], features = Crk2, split.by = 'spp', reduction = 'umap', label = T)
p2 <- FeaturePlot(S.Os[[2]], features = Crk2, split.by = 'spp', reduction = 'umap', label = T)
p3 <- FeaturePlot(S.Os[[3]], features = Crk2, split.by = 'spp', reduction = 'umap', label = T)
p4 <- FeaturePlot(S.Os[[4]], features = Crk2, split.by = 'spp', reduction = 'umap', label = T)

p <- grid.arrange(p1,p2,p3,p4, ncol = 2)

ggsave(filename="../Output/compScBdTgPb/figs/ark3_expr.pdf", plot=p,
       width = 8, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

## S phase marker expression.
toxo.markers <- read.xlsx('../Input/compScBdTgPb/gene_function/Toxo_SM_C_Markers.xlsx')


getAllFeaturePlots <- function(S.Os, my.gene){
  p1 <- FeaturePlot(S.Os[[2]], features = 'TGGT1-218520', split.by = 'spp', reduction = 'umap', label = T)
  p2 <- FeaturePlot(S.Os[[2]], features = my.gene, split.by = 'spp', reduction = 'umap', label = T)
  p3 <- FeaturePlot(S.Os[[3]], features = my.gene, split.by = 'spp', reduction = 'umap', label = T)
  p4 <- FeaturePlot(S.Os[[4]], features = my.gene, split.by = 'spp', reduction = 'umap', label = T)
  
  p <- grid.arrange(p1,p2,p3,p4, ncol = 2)
  
  return(p)
}

p <- getAllFeaturePlots(S.Os, gsub('_', '-', toxo.markers$GeneID[23]))

all.markers.sig <- left_join(all.markers.sig, prod.desc, by = 'GeneID')
all.markers.sig <- all.markers.sig %>% arrange(predicted.id, spp, desc(avg_log2FC))
write.xlsx(all.markers.sig, '../Output/compScBdTgPb/tables/up_regulated_markers_by_spp_and_phase.xlsx')


## shared PCA/UMAP
S.O.list <- SplitObject(S.Os, split.by = "spp")
S.O.list <- lapply(X = S.O.list, FUN = function(x) {
  ## Extract the count data
  
  ## extract the count data from each as.matrix(S.O.list[[1]][["RNA"]]@data)
  ## Replace genes with Bdiv orthologous when needed
  ## recreate the new Seurat object.
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = S.O.list)
reference_dataset <- ref.ind
anchors <- FindIntegrationAnchors(object.list = S.O.list, 
                                  anchor.features = features, reference = reference_dataset)
S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
S.O.integrated <- ScaleData(S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, npcs = 30, verbose = FALSE)
S.O.integrated <- RunUMAP(S.O.integrated, reduction = "pca", dims = 1:30)
S.O.integrated <- FindNeighbors(S.O.integrated, reduction = "pca", dims = 1:30)
S.O.integrated <- FindClusters(S.O.integrated, resolution = res)

S.O.integrated$phase.cond <- paste(S.O.integrated@meta.data$spp, 
                                   Idents(S.O.integrated), sep = "_")
Idents(S.O.integrated) <- "phase.cond"






## Visualize the data projected onto universal PCA
## (re) Merge the data without anchoring. This will allow the data to b visualized on
## universal PCA coordinates (so spp specific diffrences are preserved) with learned transfered
## labels.
alldata.labs <- merge(S.Os[[1]], S.Os[2:4], add.cell.ids=c("intra","extra","crk2","ark3"), merge.data = T)
alldata.labs@meta.data$spp <- factor(alldata.labs@meta.data$spp, levels = c('intra', 'extra', 'crk2', 'ark3'))
alldata.labs@meta.data$predicted.id <- factor(alldata.labs@meta.data$predicted.id, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))

alldata.labs <- prep_S.O(alldata.labs)

Idents(alldata.labs) <- 'predicted.id'
p <- DimPlot(alldata.labs, reduction = "pca", 
             #group.by = "cell", 
             split.by = 'spp',
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


alldata.labs.pcaMetadata <- getPcaMetaData(alldata.labs)

col_range <- colorRampPalette(c('white', 'yellow'))
p1 <- ggplot(alldata.labs.pcaMetadata, aes(x=UMAP_1,y=-UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = spp,
    color = spp
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('PC2') + xlab('PC1') +
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

p <- grid.arrange(p1, p2, ncol = 2)

ggsave(filename="../Output/compScBdTgPb/figs/individually_processed_universal_pca_label_transfer.pdf",
       plot=p,
       width = 12, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


library(plotly)
unique(alldata.labs.pcaMetadata$spp)
include.spps <- c('intra', 'extra', 'crk2',  'ark3')
alldata.labs.pcaMetadata.filt <- alldata.labs.pcaMetadata %>% dplyr::filter(spp %in% include.spps)
fig <- plot_ly(alldata.labs.pcaMetadata.filt, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~spp, 
               colors = colorRamp(c("green", "blue", "red", "gray66"))) 
fig <- fig %>% add_markers(size=2)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'),
                                   yaxis = list(title = 'PC_2'),
                                   zaxis = list(title = 'PC_3')))
fig




alldata <- merge(S.Os[[1]], S.Os[2:4], add.cell.ids=c("intra","extra","crk2","ark3"))
alldata@meta.data$spp <- factor(alldata@meta.data$spp, levels = c('intra', 'extra', 'crk2', 'ark3'))
alldata@meta.data$predicted.id <- factor(alldata@meta.data$predicted.id, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
alldata <- prep_S.O(alldata)

pcaMetaData.alldata <- getPcaMetaData(alldata)
pcaMetaData.alldata$predicted.id <- factor(pcaMetaData.alldata$predicted.id, 
                                   levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
pcaMetaData.alldata$phase <- factor(pcaMetaData.alldata$phase, 
                            levels = c('G1.a', 'G1.b', 'S', 'M', 'C', 'NA'))

pcaMetaData.alldata$spp <- factor(pcaMetaData.alldata$spp, 
                          levels = c('intra', 'extra', 'crk2', 'ark3'))

col_range <- colorRampPalette(c('white', 'yellow'))
p2 <- ggplot(pcaMetaData.alldata, aes(x=UMAP_1,y=-UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = predicted.id,
    color = predicted.id
  ), #color = 'blue', 
  shape=21, size = 1)+ 
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
  )
plot(p1)

p <- grid.arrange(p1, p2, ncol = 1)


## Pseudo-time
library(slingshot)
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(velociraptor)

sds <- slingshot(Embeddings(S.Os[[3]], 'umap'), 
                 clusterLabels = S.Os[[3]]@meta.data$seurat_clusters, 
                 start.clus = '1',  stretch = 2)

umap.dat <- Embeddings(S.Os[[3]], 'umap')
plot(umap.dat, pch = 20, col = 'red')
points(sds@metadata$curves$Lineage1$s, type = 'l')
points(sds@metadata$curves$Lineage2$s, type = 'l')

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sds$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
xx <- slingReducedDim(sds)
yy <- slingCurves(sds)

zz <- slingPseudotime(sds)
ww <- slingLineages(sds)



dim(xx)
sds@elementMetadata$reducedDim
sds@elementMetadata$clusterLabels

cell_colors <- cell_pal(S.Os[[3]]@meta.data$predicted.id, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())
pt <- slingPseudotime(sds)
crk.pca.umap <- getPcaMetaData(S.Os[[3]])

data.frame(sds@metadata$curves$Lineage1$s, sds@metadata$curves$Lineage1$ord, pt[!is.na(pt[,1]),1])

s.data <- data.frame(cell.ord1 = sds@metadata$curves$Lineage1$ord, cell.ord2 = sds@metadata$curves$Lineage2$ord,
                     pt = pt)
s.data$Sample <- rownames(s.data)
s.data <- s.data %>% transmute(Sample = Sample, cell.ord = cell.ord, 
                               pt = curve1, sc1 = PC_1, sc2 = PC_2)

#s.data <- left_join(s.data, S.O.combined@meta.data, by = 'Sample')


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
library(cowplot)
library(patchwork)
library(doParallel)
library(tidyverse)
library(cowplot)
library(patchwork)
library(doParallel)
library(ggVennDiagram)
library(tidytext)
library(dittoSeq)
library(ggrepel)


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

## Count files for intra and extra.
intra.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.intra.expr.csv"
extra.file.csv <- "../Input/compScBdTgPb/ToxoScRNA_MJ/RH.extra.expr.csv"

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

## generate count matrices
intra.counts <- getExpr(intra.file.csv, TGGT1_ME49)
extra.counts <- getExpr(extra.file.csv, TGGT1_ME49)


## Read in data from Lourido Alkaline stress. 
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

## 72h post-alkaline stress is the cleanest data
Meta.72 <- Meta %>% dplyr::filter(Timepoint == 72)
M.72 <- M.GT1 %>% dplyr::select(colnames(M)[colnames(M.GT1) %in% Meta.72$Sample])

M.72.WT.pH <- M.72 %>% dplyr::select(contains('WT.pH')) ## Bradyzoites
M.72.WT.D3 <- M.72 %>% dplyr::select(contains('WT.D3')) ## Tachyzoites



## Generate Seurat Objects
S.O.intra <- CreateSeuratObject(counts = intra.counts, min.cells = 10, min.features = 10)
S.O.intra$orig.ident <- 'intra'
VlnPlot(S.O.intra, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.intra, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.intra  <- subset(S.O.intra , subset = nFeature_RNA > 100 & nFeature_RNA < 1200 )

S.O.extra <- CreateSeuratObject(counts = extra.counts, min.cells = 10, min.features = 10)
S.O.extra$orig.ident <- 'extra'
VlnPlot(S.O.extra, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.extra, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.extra  <- subset(S.O.extra , subset = nFeature_RNA > 100 & nFeature_RNA < 1300 )


S.O.pru.WT <- CreateSeuratObject(counts = M.72.WT.D3, min.cells = 10, min.features = 10)
S.O.pru.WT$orig.ident <- 'pru.WT' 
VlnPlot(S.O.pru.WT, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.pru.WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.pru.WT  <- subset(S.O.pru.WT, subset = nFeature_RNA > 100 & nFeature_RNA < 1500 )

S.O.pru.WT.pH <- CreateSeuratObject(counts = M.72.WT.pH, min.cells = 10, min.features = 10)
S.O.pru.WT.pH$orig.ident <- 'pru.pH'
VlnPlot(S.O.pru.WT.pH, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.pru.WT.pH, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.pru.WT.pH  <- subset(S.O.pru.WT.pH, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 )


## Convert to list and down-sample
S.O.list <- list(intra = S.O.intra, extra = S.O.extra, pru.WT = S.O.pru.WT, pru.pH = S.O.pru.WT.pH)

S.O.list <- mclapply(S.O.list, function(S.O){
  set.seed(100)
  S.O <- subset(x = S.O, downsample = 4000)
  
}, mc.cores = num.cores)


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


saveRDS(S.Os, '../Input/compScBdTgPb/RData/S.O.intra_extra_pruWT_pruPH_lables_list.RData')

## Only intra and Extra
S.O.intra <- CreateSeuratObject(counts = intra.counts, min.cells = 10, min.features = 100)
S.O.intra$orig.ident <- 'intra'
S.O.intra$spp <- 'intra'
S.O.intra  <- subset(S.O.intra , subset = nFeature_RNA > 100 & nFeature_RNA < 1200 )
S.O.intra  <- subset(S.O.intra, downsample = 1500)

S.O.extra <- CreateSeuratObject(counts = extra.counts, min.cells = 3, min.features = 100)
S.O.extra$orig.ident <- 'extra'
S.O.extra$spp <- 'extra'
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


alldata.phase <- merge(S.Os[[1]], S.Os[2:4], add.cell.ids=c("intra","extra","pru.WT","pru.pH"))

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

#saveRDS(S.O.integrated, '../Input/compScBdTgPb/RData/S.O.intra.extra.pruWT.pruPH.integrated.rds')


## DEG analysis
DefaultAssay(S.O.integrated) <- 'RNA'

## Combine intra and pruWT for DEG analysis vs. Brady
S.O.integrated@meta.data$stage = ifelse(S.O.integrated@meta.data$spp %in% c('intra', 'pru.WT'), 'Tachy',
                                               ifelse(S.O.integrated@meta.data$spp == 'pru.pH', 'Brady', 'Extra'))
Idents(S.O.integrated) <- 'stage'
Global.Markers <- FindAllMarkers(object = S.O.integrated, only.pos = T, min.pct = 0.1) 
colnames(Global.Markers)[colnames(Global.Markers) == 'cluster'] <- 'stage'
Global.Markers$GeneID <- gsub('-', '_', Global.Markers$gene)
Global.Markers.top <- Global.Markers %>% group_by(stage) %>% slice_max(n = 1, order_by = avg_log2FC)

Global.Markers.sig <- Global.Markers %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.05) 
p <- FeaturePlot(object = S.O.integrated, features = Global.Markers.top$gene, label = T,repel = T,
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

# ggsave(filename="../Output/compScBdTgPb/figs/Tachy_Extra_Brady_total_degs.pdf", 
#        plot=p,
#        width = 6, height = 6, 
#        units = "in", # other options are "in", "cm", "mm" 
#        dpi = 300
# )


## Pairwise comparisons
## Tachy vs Brady
Idents(S.O.integrated) <- 'stage'
Brady.vs.Tachy.markers <- FindMarkers(S.O.integrated, ident.1 = "Brady", ident.2 = "Tachy", 
                                      logfc.threshold = 0, min.pct = 0.0)
Brady.vs.Tachy.markers$gene <- rownames(Brady.vs.Tachy.markers)
Brady.vs.Tachy.markers$GeneID <- gsub('-', '_', Brady.vs.Tachy.markers$gene)
Brady.vs.Tachy.markers.sig <- Brady.vs.Tachy.markers %>% dplyr::filter(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05) 
Brady.vs.Tachy.markers.sig <- left_join(Brady.vs.Tachy.markers.sig, prod.desc, by = 'GeneID') 
write.xlsx(Brady.vs.Tachy.markers.sig , '../Output/compScBdTgPb/tables/Brady_vs_Tachy_markers.xlsx')


## Tachy RH vs Extra RH 
Idents(S.O.integrated) <- 'spp'
Extra.vs.Intra.RH.markers <- FindMarkers(S.O.integrated, ident.1 = "extra", ident.2 = "intra", 
                                         logfc.threshold = 0, min.pct = 0)
Extra.vs.Intra.RH.markers$gene <- rownames(Extra.vs.Intra.RH.markers)
Extra.vs.Intra.RH.markers$GeneID <- gsub('-', '_', Extra.vs.Intra.RH.markers$gene)
Extra.vs.Intra.RH.markers.sig <- Extra.vs.Intra.RH.markers %>% dplyr::filter(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05) 
Extra.vs.Intra.RH.markers.sig <- left_join(Extra.vs.Intra.RH.markers.sig, prod.desc, by = 'GeneID') 
write.xlsx(Extra.vs.Intra.RH.markers.sig , '../Output/compScBdTgPb/tables/Extra_vs_Intra_RH_markers.xlsx')


## Shared Brady/Extra markers compared to Tachy
shared.Brady.Extra.all <- inner_join(Brady.vs.Tachy.markers, Extra.vs.Intra.RH.markers, by = 'GeneID')
shared.Brady.Extra.all$dir <- sign(shared.Brady.Extra.all$avg_log2FC.x * shared.Brady.Extra.all$avg_log2FC.y)
shared.Brady.Extra.all <- shared.Brady.Extra.all %>% arrange(desc(dir), avg_log2FC.x)
colnames(shared.Brady.Extra.all) <- gsub('\\.x', '.Brady', gsub('\\.y', '.Extra', colnames(shared.Brady.Extra.all)))
write.xlsx(shared.Brady.Extra.all, '../Output/compScBdTgPb/tables/shared_brady_extra_markers_all_genes.xlsx')

shared.Brady.Extra <- inner_join(Brady.vs.Tachy.markers.sig, Extra.vs.Intra.RH.markers.sig, by = 'GeneID')
shared.Brady.Extra$dir <- sign(shared.Brady.Extra$avg_log2FC.x * shared.Brady.Extra$avg_log2FC.y)
shared.Brady.Extra <- shared.Brady.Extra %>% arrange(desc(dir), avg_log2FC.x)
colnames(shared.Brady.Extra) <- gsub('\\.x', '.Brady', gsub('\\.y', '.Extra', colnames(shared.Brady.Extra)))
write.xlsx(shared.Brady.Extra, '../Output/compScBdTgPb/tables/shared_brady_extra_markers.xlsx')



##### Figures

S.O.intra.extra <- readRDS('../Input/compScBdTgPb/RData/S.O.intra.extra.labs_list.rds')

S.O.intra <- S.O.intra.extra$intra
S.O.intra@meta.data$spp2 <- S.O.intra@meta.data$spp
intra.pca <- getPcaMetaData(S.O.intra)

S.O.extra <- S.O.intra.extra$extra
S.O.extra@meta.data$spp2 <- S.O.extra@meta.data$spp
extra.pca <- getPcaMetaData(S.O.extra)

p1  <- ggplot(extra.pca, aes(x= -UMAP_1,y=-UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  
  theme_bw(base_size = 14) +
  #theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  #ylab('UMAP2') + xlab('UMAP1') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  # facet_grid(.~spp) + 
  # theme(
  #   strip.text.x = element_text(
  #     size = 18, color = "black", face = "bold.italic"
  #   )
  # ) + 
  theme(legend.position="none") + 
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 0),
    axis.title.y = element_text(size=18, face="bold")
  )


plot(p1)




ggsave(filename='../Output/compScBdTgPb/figs/umap_extra.pdf',
       plot=p1,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)





S.Os <- readRDS('../Input/compScBdTgPb/RData/S.O.intra_extra_pruWT_pruPH_lables_list.RData')
S.O.integrated <- readRDS('../Input/compScBdTgPb/RData/S.O.intra.extra.pruWT.pruPH.integrated.rds')



### Adding better spp names
S.O.integrated@meta.data$spp2 <- ifelse(S.O.integrated@meta.data$spp == 'intra', 'RH.intra',
                                        ifelse(S.O.integrated@meta.data$spp == 'extra', 'RH.extra',
                                               ifelse(S.O.integrated@meta.data$spp == 'pru.WT', 'Pru.intra',
                                                      'Pru.brady')))
S.O.integrated@meta.data$spp2 <- factor(S.O.integrated@meta.data$spp2, 
                                        levels = c('RH.intra', 'RH.extra', 'Pru.intra', 'Pru.brady'))
S.O.integrated$phase.spp2 <- paste(S.O.integrated@meta.data$spp2, S.O.integrated@meta.data$phase, sep = "_")
Idents(S.O.integrated) <- "phase.spp2"



integratred.pca <- getPcaMetaData(S.O.integrated)
p1  <- ggplot(integratred.pca, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  
  theme_bw(base_size = 14) +
  #theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  #ylab('UMAP2') + xlab('UMAP1') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  # facet_grid(.~spp) + 
  # theme(
  #   strip.text.x = element_text(
  #     size = 18, color = "black", face = "bold.italic"
  #   )
  # ) + 
  theme(legend.position = "None",
        #legend.position = c(0.84, 0.14),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))+ 
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 0),
    axis.title.y = element_text(size=18, face="bold")
  )


plot(p1)




ggsave(filename='../Output/compScBdTgPb/figs/umap_integrated_lourido.pdf',
       plot=p1,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



S.Os.boothroyd <- readRDS('../Input/compScBdTgPb/RData/S.O.intra_extra_pru_d0_d7_lables_list_boothroyd.RData')
S.O.integrated.boothroyd <- readRDS('../Input/compScBdTgPb/RData/S.O.intra.extra.pru_d0_d7.integrated_boothroyd.rds')

integratred.pca.boothroyd <- getPcaMetaData(S.O.integrated.boothroyd)
p1  <- ggplot(integratred.pca.boothroyd, aes(x= -UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  
  theme_bw(base_size = 14) +
  #theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  #ylab('UMAP2') + xlab('UMAP1') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  # facet_grid(.~spp) + 
  # theme(
  #   strip.text.x = element_text(
  #     size = 18, color = "black", face = "bold.italic"
  #   )
  # ) + 
  theme(legend.position = "None",
        #legend.position = c(0.84, 0.14),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))+ 
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 0),
    axis.title.y = element_text(size=18, face="bold")
  )


plot(p1)




ggsave(filename='../Output/compScBdTgPb/figs/umap_integrated_boothroyd.pdf',
       plot=p1,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## Violin plots of all known regulator
TF.info <- read.xlsx('../Input/compScBdTgPb/genes/TF_Info_Updated_kz.xlsx')
genes <- gsub("_", "-", TF.info$GeneID)
labs <- gsub("_", "-", TF.info$Ap2Name)


p <- VlnPlot(S.O.no.anchored, features = genes)

out.dir <- "../Output/compScBdTgPb/figs/expression_plots/"
p <- lapply(1:length(genes), function(i){
  pp <- p[[i]] + 
    ylab('Expr') + xlab('') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    #facet_wrap(spp~.) + 
    ggtitle(paste(labs[i], gsub('-', '_', genes[i]),sep = ':')) +
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
  
  f.name <- paste(paste(out.dir, paste(labs[i], gsub('-', '_', genes[i]),sep = '_'), sep = ''), '.pdf', sep = '')
  ggsave(filename=f.name,
         plot=pp,
         width = 5, height = 4,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
  )
  
  pp
})



## Violin plot for stage conversion players (from Sokol, Boyle review)
Stage.info <- read.xlsx('../Input/compScBdTgPb/gene_function/Brady_Stage_Conversion_Sokol_Boyle_review.xlsx')
genes <- gsub("_", "-", Stage.info$GeneID)
labs <- gsub("_", "-", Stage.info$GeneName)

DefaultAssay(S.O.integrated) <- 'RNA'
p <- VlnPlot(S.O.no.anchored, features = genes)

out.dir <- "../Output/compScBdTgPb/figs/stage_conversion_plots/"
p <- lapply(1:length(genes), function(i){
  pp <- p[[i]] + 
    ylab('Expr') + xlab('') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    #facet_wrap(spp~.) + 
    ggtitle(paste(labs[i], gsub('-', '_', genes[i]),sep = ':')) +
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
  
  f.name <- paste(paste(out.dir, paste(labs[i], gsub('-', '_', genes[i]),sep = '_'), sep = ''), '.pdf', sep = '')
  ggsave(filename=f.name,
         plot=pp,
         width = 5, height = 4,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
  )
  
  pp
})



out.dir <- '../Output/compScBdTgPb/figs/final_ap2s/'
Idents(S.O.integrated) <- 'spp2'
p <- VlnPlot(S.O.integrated, features = final.ap2.ids)
p <- lapply(1:length(final.ap2.ids), function(i){
  pp <- p[[i]] + 
    ylab('') + xlab('') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    #facet_wrap(spp~.) + 
    ggtitle(paste(final.ap2s[i], gsub('-', '_', final.ap2.ids[i]),sep = ':')) +
    theme(
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(size=12, face="bold", hjust = 0.5),
      axis.title.y = element_text(size=12, face="bold")
    ) + 
    theme(#legend.position = c(0.1, 0.25),
      legend.position = 'None',
      legend.title = element_text(colour="black", size=6, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=6, 
                                 face="bold"))
  
  f.name <- paste(paste(out.dir, paste(final.ap2s[i], gsub('-', '_', final.ap2.ids[i]),sep = '_'), sep = ''), '.pdf', sep = '')
  ggsave(filename=f.name,
         plot=pp,
         width = 4, height = 3,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
  )
  
  pp
})


## UMAP plots

pcaMataData.S.O.integrated<- getPcaMetaData(S.O.integrated)
pcaMataData.S.O.integrated$phase <- factor(pcaMataData.S.O.integrated$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))


p5  <- ggplot(pcaMataData.S.O.integrated, aes(x= UMAP_1,y=UMAP_2)) +
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
        #legend.position = c(0.88, 0.17),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=2)))


plot(p5)


ggsave(filename='../Output/compScBdTgPb/figs/projection_plots/umap_anchored_split_phase.pdf',
       plot=p5,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## markers
Global.Markers.sig <- read.xlsx('../Output/compScBdTgPb/tables/Tachy_Brady_Extra_global_markers.xlsx')

Global.Markers.stat <- Global.Markers.sig %>% group_by(stage) %>% summarise(num.deg = n())
Global.Markers.stat$stage <- factor(Global.Markers.stat$stage,
                                    levels = c('Tachy', 'Extra', 'Brady'))
# bar plot 
p <-  ggplot(Global.Markers.stat, aes(x=stage, y=num.deg)) +
  geom_bar(stat="identity", fill = "steelblue")+
  geom_text(aes(label=num.deg), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Intra (RH/Pru)", "Extra (RH)", "Brady (Pru)")) +
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

ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/Tachy_Extra_Brady_unique_markers.pdf", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


Brady.vs.Tachy.markers.sig <- read.xlsx('../Output/compScBdTgPb/tables/Brady_vs_Tachy_markers.xlsx')

Brady.vs.Tachy.markers.sig$stage <- ifelse(Brady.vs.Tachy.markers.sig$avg_log2FC > 0, 'Brady', 'Tachy')
Brady.vs.Tachy.markers.stat <- Brady.vs.Tachy.markers.sig %>% group_by(stage) %>% summarise(num.deg = n())
Brady.vs.Tachy.markers.stat$stage <- factor(Brady.vs.Tachy.markers.stat$stage,
                                            levels = c('Tachy', 'Brady'))
# bar plot 
p <-  ggplot(Brady.vs.Tachy.markers.stat, aes(x=stage, y=num.deg)) +
  geom_bar(stat="identity", fill = "steelblue")+
  geom_text(aes(label=num.deg), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Intra (RH/Pru)",  "Brady (Pru)")) +
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

ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/Tachy_Brady_total_DEGs.pdf", 
       plot=p,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

Extra.vs.Intra.RH.markers.sig <- read.xlsx('../Output/compScBdTgPb/tables/Extra_vs_Intra_RH_markers.xlsx')
Extra.vs.Intra.RH.markers.sig$stage <- ifelse(Extra.vs.Intra.RH.markers.sig$avg_log2FC > 0, 'Extra', 'Intra')
Extra.vs.Intra.RH.markers.stat <- Extra.vs.Intra.RH.markers.sig %>% group_by(stage) %>% summarise(num.deg = n())
Extra.vs.Intra.RH.markers.stat$stage <- factor(Extra.vs.Intra.RH.markers.stat$stage,
                                               levels = c('Intra', 'Extra'))
# bar plot 
p <-  ggplot(Extra.vs.Intra.RH.markers.stat, aes(x=stage, y=num.deg)) +
  geom_bar(stat="identity", fill = "steelblue")+
  geom_text(aes(label=num.deg), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Intra (RH)",  "Extra (RH)")) +
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

ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/Intra_Extra_total_DEGs.pdf", 
       plot=p,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Shared markers
shared.Brady.Extra <- read.xlsx('../Output/compScBdTgPb/tables/shared_brady_extra_markers.xlsx')
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


ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/Brady_Extra_shared_degs.pdf", 
       plot=g,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



## Up or Down Barplot
X1 <- Brady.vs.Tachy.markers.sig %>% summarise(stage = 'Brady', total = n())
X2 <- Extra.vs.Intra.RH.markers.sig %>% summarise(stage = 'Extra', total = n())
X3 <- shared.Brady.Extra %>% summarise(stage = 'Shared', total = n())
X <- rbind(X1, X2, X3)

p <-  ggplot(X, aes(x=stage, y=total, fill = stage)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=total), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Pru (pH 8.4)",  "e.c", 'shared')) +
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
        axis.ticks = element_blank()) + 
  theme(legend.position = "None",
        #legend.position = c(0.88, 0.17),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=2)))



plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/total_up_down_brady_extra_shared.pdf", 
       plot=p,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


##
shared.Brady.Extra <- read.xlsx('../Output/compScBdTgPb/tables/shared_brady_extra_markers_all_genes.xlsx')
shared.Brady.Extra <- left_join(shared.Brady.Extra, prod.desc, by = 'GeneID')

shared.Brady.Extra <- shared.Brady.Extra %>% 
  mutate(up.brady = ifelse(avg_log2FC.Brady > 1, T, F),
         up.extra = ifelse(avg_log2FC.Extra > 1, T, F))


shared.Brady.Extra <- shared.Brady.Extra %>% 
  mutate(color.class = ifelse(avg_log2FC.Brady > log2(1.7)  & avg_log2FC.Extra > log2(1.7) , 'shared', 
                              ifelse(avg_log2FC.Brady > 4 & avg_log2FC.Extra < 0.5, 'brady', 
                                     ifelse(avg_log2FC.Extra >= 2 & avg_log2FC.Brady < 0.6, 'extra','other'))))

shared.Brady.Extra$geneLable <- ''
shared.Brady.Extra <- shared.Brady.Extra %>% 
  mutate(GeneLable = ifelse(!(color.class %in%  c('other')), 'yes', 'no'))

#shared.Brady.Extra$Name <- shared.Brady.Extra$GeneID
shared.Brady.Extra$Name <- paste(shared.Brady.Extra$GeneID, shared.Brady.Extra$ProductDescription, sep = '\n')

shared.Brady.Extra$Name <- gsub('kinase', 'kin.', gsub('serine/threonine prot. phosphatase', 'ser./thre. phosph. ', gsub('Rhoptry kinase fam. prot.', 'Rhop. kinase', gsub('tryptophanyl-tRNA synthetase \\(TrpRS2\\)', 'TrpRS2', 
                                                                                                                                                                           gsub('putative cell-cycle-associated protein kinase PRP4', 'PRP4', 
                                                                                                                                                                                gsub('hop. kinase fam. prot.', 'Rhop. kinase ', 
                                                                                                                                                                                     gsub('AP2 domain-containing protein', 'AP2 dom.', 
                                                                                                                                                                                          gsub('bradyzoite rhoptry protein ', '', 
                                                                                                                                                                                               gsub('threonine specific protein phosphatase', 'threonine prot. phosphatase', 
                                                                                                                                                                                                    gsub('cAMP-dependent protein kinase', 'cAMP-dep kinase', 
                                                                                                                                                                                                         gsub('universal', 'univ.', gsub('phosphatidylinositol 3- and 4-kinase', 'kinase', 
                                                                                                                                                                                                                                         gsub('putative cell-cycle-associated protein kinase GSK', 'kinase GSK', 
                                                                                                                                                                                                                                              gsub('bradyzoite antigen ', '', 
                                                                                                                                                                                                                                                   gsub('putative rhoptry kinase fam. prot. ROP34', 'ROP34', 
                                                                                                                                                                                                                                                        gsub('protease inhibitor PI2', 'protease inh. PI2', 
                                                                                                                                                                                                                                                             gsub('subtilisin SUB1', 'SUB1', 
                                                                                                                                                                                                                                                                  gsub('lactate dehydrogenase ', '', 
                                                                                                                                                                                                                                                                       gsub('tetratricopeptide repeat-containing protein', 'tetratricopeptide', 
                                                                                                                                                                                                                                                                            gsub('putative microneme protein', 'microneme', 
                                                                                                                                                                                                                                                                                 gsub('AP2 domain transcription factor ', '', 
                                                                                                                                                                                                                                                                                      gsub('putative myosin heavy chain', 'GRA46',
                                                                                                                                                                                                                                                                                           gsub('zinc finger \\(CCCH type\\) motif-containing protein', 'Zn fing.', 
                                                                                                                                                                                                                                                                                                gsub('transporter', 'trans.', gsub('SAG-related sequence ', '', 
                                                                                                                                                                                                                                                                                                                                   gsub('family protein', 'fam. prot.', 
                                                                                                                                                                                                                                                                                                                                        gsub('OTU family cysteine protease', 'cysteine protease', 
                                                                                                                                                                                                                                                                                                                                             gsub('hypothetical protein', 'hypo.',
                                                                                                                                                                                                                                                                                                                                                  gsub('tubulin/FtsZ family, GTPase domain-containing protein', 'tub epsilon',
                                                                                                                                                                                                                                                                                                                                                       shared.Brady.Extra$Name)))))))))))))))))))))))))))))

shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_311100')] <- 'TGGT1_311100 \n BFD2'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_215210')] <- 'TGGT1_215210 \n F-Box'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_239748')] <- 'TGGT1_239748 \n GPI-GlcNAc transf. \n (fitness!)'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_215970')] <- 'TGGT1_215970 \n hypo. \n excreted'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_247530')] <- 'TGGT1_247530 \n hypo. 1xTM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_203300')] <- 'TGGT1_203300 \n hypo. 2xTM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_217530')] <- 'TGGT1_217530 \n GRA63'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_265320')] <- 'TGGT1_265320 \n splicing factor CWC2'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_311230')] <- 'TGGT1_311230 \n mucin? TM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_273320')] <- 'TGGT1_273320 \n hypo. SP & TM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_321530')] <- 'TGGT1_321530 \n Catheptsin L'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_276170')] <- 'TGGT1_276170 \n PI3,4K'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_312330')] <- 'TGGT1_312330 \n RNA-PolII-B1'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_207160')] <- 'TGGT1_207160 \n SRS49D'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_268860')] <- 'TGGT1_268860 \n ENO1'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_248990')] <- 'TGGT1_248990 \n put. prenylcysteine oxidase'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_272370')] <- 'TGGT1_272370 \n hypo. Sx TM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_223550')] <- 'TGGT1_223550 \n hypo.I 4xTM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_240090')] <- 'TGGT1_240090 \n ROP34'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_208730')] <- 'TGGT1_208730 \n MIC (putative)'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_225540')] <- 'TGGT1_225540 \n hypo. (TolA?)'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_216140')] <- 'TGGT1_216140 \n TPR domain \n (pub)'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_269690')] <- 'TGGT1_269690 \n hypo. SP'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_253330')] <- 'TGGT1_253330 \n BPK1 (ROP kin)'












shared.Brady.Extra$Name[shared.Brady.Extra$GeneLable == 'yes']
p <- ggplot(shared.Brady.Extra) + geom_point(aes(avg_log2FC.Brady, avg_log2FC.Extra,color=color.class))+ 
  #xlim(-4, 4) + ylim(-3, 3) + 
  labs(color = "") + 
  geom_text_repel(aes(avg_log2FC.Brady, avg_log2FC.Extra), size = 1.7, fontface = "bold",
                  label = ifelse(shared.Brady.Extra$GeneLable == 'yes', as.character(shared.Brady.Extra$Name),""), 
                  box.padding = unit(0.6, "lines"),
                  max.overlaps = 300,
                  segment.angle = 180,
                  nudge_x = 0.5, 
                  nudge_y = 0.5,
                  hjust=0,
                  #nudge_x=0.25, 
                  segment.size = 0.25) + 
  geom_hline(yintercept = log2(1.7), color = 'black', linetype=2) + 
  geom_vline(xintercept = log2(1.7), color = 'black', linetype=2) + 
  xlab("log2 FC Brady over Intra") + ylab("log2 FC Extra over Intra") +
  theme(legend.title=element_blank(),text = element_text(size=20))+ 
  scale_color_manual(values = c("shared" = "red", "brady" = "blue", 'extra' = 'green', 'other' = 'gray'))+
  theme_bw()+
  theme(
    #plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    legend.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.45, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    strip.text = element_text(color = "black", size = 18, face = "bold"),
    strip.background = element_rect(fill = "white"))+
  theme(legend.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.position = c(0.08, 0.25))

p

ggsave(filename="../Output/compScBdTgPb/figs/brady_over_intra_vs_extra_over_intra_fc1_7.pdf", 
       plot=p,
       width = 9, height = 9, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



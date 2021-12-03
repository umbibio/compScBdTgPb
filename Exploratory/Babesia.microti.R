library(Seurat)
library(dplyr)
library(openxlsx)
library(tidyverse)
# Load the dataset
Bm <- Read10X(data.dir = "../Input/Bmicroti/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
Bm <- CreateSeuratObject(counts = Bm, project = "NAME", min.cells = 3, min.features = 10)
Bm

#calculate the number of drop outs 
mat<-Bm[["RNA"]]@counts
sum(mat == 0) / (dim(mat)[1]*dim(mat)[2] )
Bm[["percent.mt"]] <- PercentageFeatureSet(object = Bm, pattern = "^BmR1-mt")
VlnPlot(Bm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Bm <- subset(x = Bm, subset =  percent.mt < 1)
#Bm

plot1 <- FeatureScatter(object = Bm, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot1

plot2 <- FeatureScatter(object = Bm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

feature<- Bm[["nFeature_RNA"]]
summary(feature)
counts <- Bm[["nCount_RNA"]]
summary(counts)
boxplot(feature, main="Distribution of Feature Counts")
boxplot(counts, main="Distribution of Raw Counts")

VlnPlot(Bm, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

unlistFeat <-unlist(feature)
numeric <- as.numeric(unlistFeat)
sd(numeric, na.rm = FALSE)

Bm <- subset(x = Bm, subset = nFeature_RNA > 50 & nFeature_RNA < 300 &  percent.mt < 3)
Bm
Bm <- NormalizeData(object = Bm)

Bm <- FindVariableFeatures(object = Bm, selection.method = 'vst', nfeatures = 2000)

top10 <- head(x = VariableFeatures(object = Bm), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = Bm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot1
plot2

top10
#nowscale the genes
all.genes <- rownames(x = Bm)
Bm <- ScaleData(object = Bm , features = all.genes)

#nwo run PCA dimension reduction
Bm <- RunPCA(object = Bm, features = VariableFeatures(object = Bm))


#find significant PCA by xonstructing elbow plot. 
ElbowPlot(Bm)

#visualize PCA 


print(x = Bm[['pca']], dims = 1:7, nfeatures = 10)
VizDimLoadings(object = Bm, dims = 1:7, reduction = 'pca')
DimPlot(object = Bm, reduction = 'pca')

#show the top 5 genes in the most significant PCs
print(Bm[["pca"]], dims = 1:7, nfeatures = 5)

#run jackstraw plot to see how many PCs are important 
Bm <- JackStraw(Bm, num.replicate = 100)
Bm <- ScoreJackStraw(Bm, dims = 1:20)
JackStrawPlot(Bm, dims = 1:20)

# I am going to run clustering on the top 4 PCs 
#clustering 
Bm <- FindNeighbors(Bm, dims = 1:7)
Bm <- FindClusters(Bm, resolution = 0.2)


Bm <- RunUMAP(object = Bm, dims = 1:7)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(object = Bm, reduction = 'umap')

###


##Get Sexual Stage gene
Sextual.Stage <- read.xlsx('../Input/SextualStageMarkers/Babesia_Sexual_Stage_Markers_KZ.xlsx', sheet = 3)



prod.desc <- read.csv('../Input/BdivCellCycle/ProductDescription/BDvi_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

#Idents(all.samples.integrated) <- "phase.cond"
p1 <- DimPlot(Bm, reduction = "pca", dims = c(1,2),
              #group.by = "cells", 
              #split.by = 'spp',
              pt.size = 1,
              #shape.by='spp',
              label = TRUE, label.size = 6) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p1)


ggsave(filename="../Output/scClockFigs/b_microti_pca_filter_mitochondria.pdf",
       plot=p1,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600)
       


### Check sexual stage marker expression

sextual.marker <- gsub(" ", "", gsub('_', '-', Sextual.Stage$Bm))
p3 <- FeaturePlot(object = Bm, 
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = sextual.marker,
                  cols = c("lightgrey", "red"), reduction = "umap") 

plot(p3)

ggsave(filename="../Output/scClockFigs/b_microti_sextual_markers.pdf",
       plot=p3,
       width = 12, height = 12,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)


p4 <- FeaturePlot(object = Bm, 
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = sextual.marker[c(1,2)], blend = T, 
                  cols = c("lightgrey", "red", 'blue'), reduction = "umap") 

ggsave(filename="../Output/scClockFigs/b_microti_sextual_markers_blend_1_2.pdf",
       plot=p4[[3]],
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)



Bm@meta.data$Sample <- rownames(Bm@meta.data)
sexual.expressing.ind <- which(Bm[["RNA"]]@data[grep(sextual.marker[2], rownames(Bm[["RNA"]]@data)),] > 0)
length(sexual.expressing.ind)
sextual.expressing <- colnames(Bm[["RNA"]]@data)[sexual.expressing.ind]
Bm@meta.data$sextual <- ifelse(Bm@meta.data$Sample %in% sextual.expressing, 1, 0)
Idents(Bm) <- 'sextual'

Sextual.markers <- FindAllMarkers(object = Bm, only.pos = TRUE, min.pct = 0)

Sextual.markers$GeneID <- gsub('-', '_', Sextual.markers$gene)
Sextual.markers.sig <- Sextual.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05 & cluster == 1)
sextual.marker %in% Sextual.markers.sig$gene

write.xlsx(PKG.markers.sig, '../Output/scClockOut/PKG_expressing_cells_markers.xlsx')
FeaturePlot(object = S.O.bd.filt, 
            features = PKG.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")






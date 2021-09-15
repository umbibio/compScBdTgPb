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
#library(sctransform)


source('./util_funcs.R')


### T. gondii
input.dir.tg <- '../Input/scRNAseqToxoAtlas/kz/'

rh384.expr.file <- 'rh384_expression_filtered.csv'
rh384.obs.file <- 'rh384_obs_filtered.csv'


## Reading T. gondii data
rh384.expr <- read.csv(paste(input.dir.tg, rh384.expr.file, sep = ''))
rh384.obs <- read.csv(paste(input.dir.tg, rh384.obs.file, sep = ''))

processCounts <- function(expr){
  cols <- expr[,1]
  expr <- t(expr[,-1])
  colnames(expr) <- cols
  return(expr)
}

rh384.expr <- processCounts(rh384.expr)
rh384.obs$spp <- 'Tg'

rh384.obs$phase <- gsub(' ', '.', gsub("\"", "",  as.character(rh384.obs$sub_cell_cycle)))
rh384.obs$cells <- paste('Tg', gsub("\"", "",  as.character(rh384.obs$cell_cycle)), sep = '.')
rh384.obs$NAME <- paste(rh384.obs$cells, rh384.obs$index, sep = '_')
tg.pheno <- rh384.obs %>% transmute(Sample = index, spp = spp, cells = cells, Name = NAME, phase = phase)

# Set initial Seurat clusters
S.O.tg <- CreateSeuratObject(counts = rh384.expr, min.cells = 10, min.features = 100)


#VlnPlot(S.O.tg, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O.tg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

S.O.tg <- subset(S.O.tg, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
S.O.tg <- prep_S.O(S.O.tg)
tg.pheno <- tg.pheno %>% dplyr::filter(Sample %in% colnames(S.O.tg))
rownames(tg.pheno) <- tg.pheno$Sample
S.O.tg <- AddMetaData(S.O.tg, tg.pheno)


saveRDS(S.O.tg, '../Input/compScBdTgPb/RData/S.O.tg.RData')

# ## Differential gene expression
Idents(S.O.tg) <- 'phase'
Tg.markers <- FindAllMarkers(object = S.O.tg, only.pos = TRUE, min.pct = 0)

Tg.markers$GeneID <- gsub('-', '_', Tg.markers$gene)
Tg.markers.top <- Tg.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = S.O.tg, 
            features = Tg.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

Tg.markers.sig <- Tg.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)

saveRDS(Tg.markers.sig, '../Input/compScBdTgPb/RData/TG.markers.sig.RData')
ss <- Tg.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
ss$cluster <- factor(ss$cluster, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
p <- ggplot(data=ss, aes(x=cluster, y=num.DEG)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="white", size=3.5)+
  theme_minimal()

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/tg_deg_numbers.pdf",
       plot=p,
       width = 8, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

print(ss)


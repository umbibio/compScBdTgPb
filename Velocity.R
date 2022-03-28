library(velocyto.R)
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library(openxlsx)

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

## loom files generated using velocyto
## On Markov:
## source .bashrc
## conda activate RNAvel
## velocyto run10x /H3/scRNAseqBabesCS/count/bd7dOUT_count/ /H3/scRNAseqBabesCS/Genome/PiroplasmaDB-55_Bdivergens1802A.gtf
## conda deactivate



input.dir.bdiv.loom <- "../Input/compScBdTgPb/loom/"
loom.files <- list.files(input.dir.bdiv.loom)
num.total.files <- length(loom.files)

loom.list <- list()

for(i in 1:num.total.files){
  cat(paste('processing file', loom.files[i]))
  cat('\n')
  ldat <- read.loom.matrices(paste(input.dir.bdiv.loom, loom.files[i], sep = ''))
  loom.list <- c(loom.list, list(ldat))
}

names(loom.list) <- c('intra','extra', 'ark3','crk2')

loom.list <- lapply(loom.list, function(ldat){
  ldat <- lapply(ldat,function(x) {
    colnames(x) <-  gsub("x","\\.1",gsub(".*:","",colnames(x)))
    x
  })
})


#S.Os <- readRDS('../Input/compScBdTgPb/RData/S.O.intra_extra_crk2_ark3_lables_not_anchored.RData')

prod.desc  <- read.xlsx('../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## Data from Boothroyd for label transfer
S.O.tg <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.RData')

spps <- names(loom.list)

prep.Loom <- function(ldat){
  
  ## Go through assays (spliced, unspliced ambiguous) and change Gene IDs to TGGT1
  ldat <- lapply(ldat, function(x){
    genes <- rownames(x)
    ind <- which(genes %in% TGGT1_ME49$TGME49)
    counts <- x[ind, ]
    genes <- genes[ind]
    genes <- TGGT1_ME49$TGGT1[match(genes, TGGT1_ME49$TGME49)]
    rownames(counts) <- genes
    counts
  })  
  
  bm <- as.Seurat(x = ldat)
  
  #VlnPlot(bm, features = c("nFeature_spliced", "nCount_spliced"), ncol = 2)
  #FeatureScatter(bm, feature1 = "nCount_spliced", feature2 = "nFeature_spliced")
  bm <- subset(bm, subset = nFeature_spliced > 100 & nFeature_spliced < 1200 )
  
  bm <- NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
  bm <- FindVariableFeatures(bm, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(bm)
  bm <- ScaleData(bm, features = all.genes)
  #bm <- SCTransform(object = bm, assay = "spliced")
  bm <- RunPCA(object = bm, verbose = FALSE)
  bm <- FindNeighbors(object = bm, dims = 1:20)
  bm <- FindClusters(object = bm)
  bm <- RunUMAP(object = bm, dims = 1:20)
  
  ## Transfer labels from Boothroyd
  anchors <- FindTransferAnchors(reference = S.O.tg, query = bm, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.tg@meta.data$phase,dims = 1:30)
  predictions$phase <- predictions$predicted.id
  bm <- AddMetaData(object = bm, metadata = predictions)
  bm@meta.data$phase <- factor(as.character(bm@meta.data$phase), levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
  
  bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
  return(bm)
}


bms <- lapply(loom.list, function(ldat){
  bm <- prep.Loom(ldat)
  Idents(bm) <- 'phase'
  
  return(bm)
})

saveRDS(bms, '../Input/compScBdTgPb/RData/rna_vel_loom_toxo_labs.RData')

#ident.colors <- (scales::hue_pal())(n = length(c('G1.a', 'G1.b', 'S', 'M', 'C')))
ident.colors <- c("#af2225",'#ee7400', '#caae05', '#6f893e', '#b437f1')
names(x = ident.colors) <- c('G1.a', 'G1.b', 'S', 'M', 'C')

cell.colors <- lapply(bms, function(bm){
  cell.color <- ident.colors[Idents(object = bm)]
  names(x = cell.color) <- colnames(x = bm)
  cell.color
})




pdf("../Output/compScBdTgPb/figs/intra_umap_vel.pdf",
    width = 6, height = 6
) 


par(mar = c(5, 5, 4, 4) + 0.1)
show.velocity.on.embedding.cor(emb = Embeddings(object = bms$intra, reduction = "umap"), 
                               vel = Tool(object = bms$intra, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors$intra, alpha = 0.5), 
                               cex = 1.8, arrow.scale = 2, show.grid.flow = TRUE, xlab = 'UMAP1', ylab = 'UMAP2',
                               cex.lab = 1.8, cex.axis = 1.2,
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

dev.off()




pdf("../Output/compScBdTgPb/figs/extra_umap_vel.pdf",
    width = 6, height = 6
) 


par(mar = c(5, 5, 4, 4) + 0.1)
show.velocity.on.embedding.cor(emb = Embeddings(object = bms$extra, reduction = "umap"), 
                               vel = Tool(object = bms$extra, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors$extra, alpha = 0.5), 
                               cex = 1.8, arrow.scale = 2, show.grid.flow = TRUE, xlab = 'UMAP1', ylab = 'UMAP2',
                               cex.lab = 1.8, cex.axis = 1.2,
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

dev.off()



pdf("../Output/compScBdTgPb/figs/crk2_umap_vel.pdf",
    width = 6, height = 6
) 


par(mar = c(5, 5, 4, 4) + 0.1)
show.velocity.on.embedding.cor(emb = Embeddings(object = bms$crk2, reduction = "umap"), 
                               vel = Tool(object = bms$crk2, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors$crk2, alpha = 0.5), 
                               cex = 1.8, arrow.scale = 2, show.grid.flow = TRUE, xlab = 'UMAP1', ylab = 'UMAP2',
                               cex.lab = 1.8, cex.axis = 1.2,
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

dev.off()


pdf("../Output/compScBdTgPb/figs/ark3_umap_vel.pdf",
    width = 6, height = 6
) 

par(mar = c(5, 5, 4, 4) + 0.1)
show.velocity.on.embedding.cor(emb = Embeddings(object = bms$ark3, reduction = "umap"), 
                               vel = Tool(object = bms$ark3, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors$ark3, alpha = 0.5), 
                               cex = 1.8, arrow.scale = 2, show.grid.flow = TRUE, xlab = 'UMAP1', ylab = 'UMAP2',
                               cex.lab = 1.8, cex.axis = 1.2,
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

dev.off()




##### Old stuff
library(velocyto.R)
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
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

## loom files generated using velocyto
## On Markov:
## conda activate RNAvel
## velocyto run10x /H3/scRNAseqToxo/NextSeq/count_ME49/6-hr-extracellular_RH-YFP2_ME49_count/ /H3/scRNAseqToxo/NextSeq/Genome/ME49_Genome/ToxoDB-54_TgondiiME49.gtf

ldat <- read.loom.matrices("../Input/compScBdTgPb/loom/27-30-33-hpi_with_10M_RH-YFP2_ME49_count.loom")
ldat <- read.loom.matrices("../Input/compScBdTgPb/loom/deltaCrk2_2-hr-IAA-treated_ME49_count.loom")
ldat <- read.loom.matrices("../Input/compScBdTgPb/loom/deltaArk3_6-hr-ATc-treated_ME49_count.loom")
ldat <- read.loom.matrices("../Input/compScBdTgPb/loom/6-hr-extracellular_RH-YFP2_ME49_count.loom")
str(ldat)
#colnames(ldat[[1]])

ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("x","\\.1",gsub(".*:","",colnames(x)))
  x
})


S.O <- readRDS('../Input/compScBdTgPb/RData/S.O.RH.crk2.RData')

pca.dat <- getPcaMetaData(S.O)




bm <- as.Seurat(x = ldat)
bm <- SCTransform(object = bm, assay = "spliced")
bm <- RunPCA(object = bm, verbose = FALSE)
bm <- FindNeighbors(object = bm, dims = 1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm, dims = 1:20)
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), 
                               vel = Tool(object = bm, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


emb <- pca.dat %>% select(PC_1, PC_2)
rownames(emb) <- pca.dat$Sample

cell.colors <- pca.dat %>% select(predicted.id)

cell.colors <- cell.colors %>% 
  mutate(color = case_when(predicted.id == "G1.a" ~ "firebrick",
                           predicted.id == "G1.b" ~ "darkorange2", 
                           predicted.id == 'S' ~ 'gold3', 
                           predicted.id == 'M' ~ 'darkolivegreen4', 
                           predicted.id == 'C' ~ 'darkorchid2'))
rownames(cell.colors) <- pca.dat$Sample
cc <- cell.colors$color
names(cc) <- rownames(cell.colors)


cell.colors <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/cell.colors.rds"))
emb <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/embedding.rds"))


hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning
# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,cc,min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat,cc,min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat,cc,min.max.cluster.average = 0.8)
# look at the resulting gene set
length(intersect(rownames(emat),rownames(nmat)))


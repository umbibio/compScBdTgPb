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

source('./util_funcs.R')

## This script reads the expr.csv files for all Babesia species
## It generates and saves the Seurat object. 

input.dir <- "../Input/compScBdTgPb/ToxoScRNA_MJ/"
count.files <- list.files(input.dir)
num.total.files <- length(count.files)
  
file.info <- data.frame(Species = rep("RH", num.total.files))

## These were determined with visual inspection
RH.intra.cutoff <- c(200, 1500)
RH.extra.cutoff <- c(200, 1500)
RH.ark3.cutoff  <- c(200, 1500)
RH.crk2.cutoff  <- c(200, 1500)

processCountToxo <- function(input.dir, filename, down.sample = T){
  file.csv <- paste(filename, ".expr.csv", sep = "")
  
  print(file.csv)
  
  file.counts <- read.csv(paste(input.dir, file.csv, sep = '/'))
  genes <- file.counts$X
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  
  spp <- filename
  print(spp)
  if(grepl('INTRA', toupper(filename))){
    cutoff <- RH.intra.cutoff
  }else if(grepl('EXTRA', toupper(filename))){
    cutoff <- RH.extra.cutoff
   }else if(grepl('ARK3', toupper(filename))){
    cutoff <- RH.ark3.cutoff
  }else if(grepl('CRK2', toupper(filename))){
    cutoff <- RH.crk2.cutoff
  }
  
  print(cutoff)

  S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
  S.O <- subset(S.O, subset = nFeature_RNA > cutoff[1] & nFeature_RNA < cutoff[2])
  
  
  
  #VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.01, 0.9))
  #print(cutoffs)
  #S.O <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )
  
  
  S.O <- prep_S.O(S.O)
  
  
  cluster <- as.character(S.O$seurat_clusters)
  
  pheno <- data.frame(Sample = names(S.O$orig.ident))
  pheno$spp <- spp
  pheno$cluster <- cluster
  pheno$cell <- paste(pheno$spp, pheno$cluster , sep = "")
  pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = "_")
  
  ind <- match(pheno$Sample,  colnames(expr))
  colnames(expr) <- pheno$NAME[ind]
  
  if(down.sample){
    set.seed(100)
    S.O <- subset(x = S.O, downsample = 800)
  }
  
  S.O.obj <- paste(paste("S.O", spp, sep = "."), "RData", sep = ".")
  S.O.dir <- paste('../Input/compScBdTgPb/RData/' , S.O.obj, sep = "")
  saveRDS(S.O, S.O.dir)
  
  L <- list(pheno = pheno, S.O = S.O)
  return(L)
  
  return(L)
}

S.O.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- gsub('.expr.csv', '', count.files[i])
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  L <- processCountToxo(input.dir, file.info$filename[i])
  S.O.list <- c(S.O.list, list(L))
}

names(S.O.list) <- c('RH_ark3', 'RH_crk2','RH_extra', 'RH_intra')
saveRDS(S.O.list, '../Input/compScBdTgPb//RData/S.O.Toxo.list.Rdata')


### Working with orthologous genes only
orthologs <- read.xlsx("../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx")
#orthologs <- read.xlsx("../Input/compScBabesia/Orthologs/Bdiv_Bbig_Bbov_orth.xlsx")


num.objs <- length(S.O.list)

phenos <- lapply(S.O.list, `[[`, 1)
S.Os <-  lapply(S.O.list, `[[`, 2)

spps <- unlist(lapply(phenos, function(x) x$spp[1]))


ortho.genes.GT1.id <- lapply(1:num.objs, function(i){
  ind.orth <- which(gsub('-', '_', rownames(S.Os[[i]]@assays$RNA@data)) %in% orthologs$TGME49)
  genes.with.orth <- gsub('-', '_', rownames(S.Os[[i]]@assays$RNA@data))[ind.orth]
  genes.GT1.id <- orthologs$TGGT1[match(genes.with.orth,orthologs$TGME49)]
} )

common.genes <- Reduce(intersect, ortho.genes.GT1.id)
orthologs.common <- orthologs[which(orthologs$TGGT1 %in% common.genes), ]


processCountToxoOrth <- function(input.dir, filename, orthologs.common, down.sample = T){
  file.csv <- paste(filename, ".expr.csv", sep = "")
  
  file.counts <- read.csv(paste(input.dir, file.csv, sep = '/'))
  genes <- file.counts$X
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  expr.mat <- as.matrix(S.Os[[i]]@assays$RNA@data)
  
  spp <- file.info$filename[i]
  
  expr <- expr[which(gsub('-', '_',rownames(expr)) %in% orthologs.common$TGME49),]
  rownames(expr) <- orthologs.common$TGGT1[match(gsub('-', '_',rownames(expr)), orthologs.common$TGME49)]
  
  spp <- filename
  print(spp)
  if(grepl('INTRA', toupper(filename))){
    cutoff <- RH.intra.cutoff
  }else if(grepl('EXTRA', toupper(filename))){
    cutoff <- RH.extra.cutoff
  }else if(grepl('ARK3', toupper(filename))){
    cutoff <- RH.ark3.cutoff
  }else if(grepl('CRK2', toupper(filename))){
    cutoff <- RH.crk2.cutoff
  }
  
  print(cutoff)
  
  S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
  S.O <- subset(S.O, subset = nFeature_RNA > cutoff[1] & nFeature_RNA < cutoff[2])
  
  
  #FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.01, 0.9))
  #print(cutoffs)
  #S.O <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )
  
  
  S.O <- prep_S.O(S.O)
  
  
  cluster <- as.character(S.O$seurat_clusters)
  
  pheno <- data.frame(Sample = names(S.O$orig.ident))
  pheno$spp <- spp
  pheno$cluster <- cluster
  pheno$cell <- paste(pheno$spp, pheno$cluster , sep = "")
  pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = "_")
  
  ind <- match(pheno$Sample,  colnames(expr))
  colnames(expr) <- pheno$NAME[ind]
  
  if(down.sample){
    set.seed(100)
    S.O <- subset(x = S.O, downsample = 800)
  }
  
  S.O.obj <- paste(paste("S.O", spp, "ortholog", sep = "."), "RData", sep = ".")
  S.O.dir <- paste('../Input/compScBabesia/RData/' , S.O.obj, sep = "")
  saveRDS(S.O, S.O.dir)
  
  L <- list(pheno = pheno, S.O = S.O)
  return(L)
}

S.O.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- gsub('.expr.csv', '', count.files[i])
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  L <- processCountToxoOrth(input.dir, file.info$filename[i], orthologs.common)
  S.O.list <- c(S.O.list, list(L))
}

names(S.O.list) <- c('RH_ark3', 'RH_crk2','RH_extra', 'RH_intra')
saveRDS(S.O.list, '../Input/compScBdTgPb//RData/S.O.Toxo.list.ortholog.Rdata')

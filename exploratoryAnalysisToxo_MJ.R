library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(sctransform)

source('./util_funcs.R')

S.O.list <- readRDS('../Input/compScBdTgPb/RData/S.O.Toxo.list.ortholog.Rdata')
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

phenos <- lapply(S.O.list, `[[`, 1)
S.Os <-  lapply(S.O.list, `[[`, 2)

spps <- unlist(lapply(phenos, function(x) x$spp[1]))
print(spps)

## Integrate Samples
ref.ind <- 4
data.ind <- c(1,2,3,4) 

all.samples.integrated <- processeMergedS.O(S.O.list, file.info, data.ind, ref.in = NA, res = 0.2, SC = F)

all.samples.integrated@meta.data$spp <- factor(all.samples.integrated@meta.data$spp, 
                                               levels = unique(all.samples.integrated@meta.data$spp))

saveRDS(all.samples.integrated, '../Input/compScBabesia/RData/all.samples.integrated.RData')

Idents(all.samples.integrated) <- "cluster"

p1 <- DimPlot(all.samples.integrated, reduction = "umap", group.by = "spp")
p2 <- DimPlot(all.samples.integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DefaultAssay(all.samples.integrated) <- 'integrated'
p <- DimPlot(all.samples.integrated, reduction = "umap", 
             #group.by = "spp", 
             split.by = 'spp',
             pt.size = 1,
             #shape.by='spp',
             label = TRUE, label.size = 4) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)

ggsave(filename="../Output/compScBabsSpecies/figs/toxo_new_integrated_pca.png",
       plot=p,
       width = 12, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)






DefaultAssay(all.samples.integrated) <- "RNA"
Idents(all.samples.integrated) <- "spp"

markers.int <- FindAllMarkers(object = all.samples.integrated, only.pos = TRUE, min.pct = 0.25) 
colnames(markers.int)[colnames(markers.int) == 'cluster'] <- 'spp'
markers.int$GeneID <- gsub('-', '_', markers.int$gene)
markers.int.sig <- markers.int %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
markers.int.top <- markers.int.sig %>% group_by(spp) %>% slice_max(n = 1, order_by = avg_log2FC)


#### Here we need to generate a bar plot of spp specific global markers total numbers.
markers.int.sig.stats <- markers.int.sig %>% group_by(spp) %>% summarise(num.deg = n())
markers.int.sig.stats$spp <- factor(markers.int.sig.stats$spp, levels = c('RH.intra', 'RH.extra', 'RH.ark3', 'RH.crk2'))
p <- ggplot(data=markers.int.sig.stats, aes(x=spp, y=num.deg)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.deg), vjust=1.6, color="white", size=3.5)+
  theme_minimal()


plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/tg_new_deg_numbers.pdf",
       plot=p,
       width = 8, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


prod.desc <- read.xlsx('../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx')
ME49.GT1 <- read.xlsx('../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, ME49.GT1, by = c('GeneID' = 'TGGT1'))
markers.int.sig <- left_join(markers.int.sig, prod.desc, by = c('GeneID' = 'TGME49')) %>% arrange(spp, desc(avg_log2FC))


write.xlsx(markers.int.sig, '../Output/compScBdTgPb/tables/RH_integrated_data_uniq_glob_marker.xlsx')




DefaultAssay(all.samples.integrated) <- "RNA"
Idents(all.samples.integrated) <- 'cluster'
nk.markers <- FindConservedMarkers(all.samples.integrated, ident.1 = 0, grouping.var = "spp", verbose = T)
head(nk.markers)




## Split the integrated S.Os
DefaultAssay(all.samples.integrated) <- "RNA" ## switch back to original RNA assay

S.O.integrated.list <- SplitObject(all.samples.integrated, split.by = 'spp')

## Fit a pseudo-time to each using and align with Bdiv bulk data.
bd.tc.logCPM <- readRDS('../Input/compScBabesia/RData/bd_sync_tc_logCPM.RData')

plotUpdatedPstimeS.O <- function(L){
  Idents(L$S.O.bd.update) <- 'adj.time.idx'
  p <- DimPlot(L$S.O.bd.update, reduction = "pca", 
               #group.by = "cells", 
               #split.by = 'spp',
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
  
  return(p)
  
}

getMaximaLocs <- function(t, y){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) > 1){
    maxima.ind = {}
    for(i in inc.ind){
      maxima.ind = c(maxima.ind, sum(locs$lengths[1:i]))
    }
    ## Interpolate a point between the location where derivative changes sign
    maxima = (spline.fit$x[maxima.ind] + spline.fit$x[(maxima.ind + 1)]) / 2
    maxima = maxima[!is.na(maxima)]
    ## Get the maximum values
    maxval = predict(spline.fit, maxima)
    
    ## Get the outliers
    #maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = 0.8))
    
    ## Peaks for entities of interest
    entity.x = maxval$x
    entity.y = maxval$y
    
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  L <- list(entity.x = entity.x, entity.y = entity.y)
  
  return(L)
}

### Fit pseudo-time, align with B div bulk, and manually adjust based on distribution of lag times
print(spps)

L1 <- list()
L2 <- list()

## BBIG
L1$bbig <- fitPseudoTime(S.O.integrated.list[[1]], reverset.time = T)
L2$bbig <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[1]], L1$bbig$cell.cycle.genes.df, 
                                           L1$bbig$sds.data, bd.tc.logCPM)

plot(L2$bbig$den)
print(L2$bbig$lag.time)
p <- plotUpdatedPstimeS.O(L2$bbig)
plot(p)


## Adjusting the start
getMaximaLocs(L2$bbig$den$x, L2$bbig$den$y) 

L2$bbig <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[1]], L1$bbig$cell.cycle.genes.df, 
                                           L1$bbig$sds.data, bd.tc.logCPM, lag.time = 30)

p <- plotUpdatedPstimeS.O(L2$bbig)
plot(p)


## BBOV
L1$bbov <- fitPseudoTime(S.O.integrated.list[[2]], reverset.time = T)
L2$bbov <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[2]], L1$bbov$cell.cycle.genes.df, 
                                           L1$bbov$sds.data, bd.tc.logCPM)

plot(L2$bbov$den)
print(L2$bbov$lag.time)
p <- plotUpdatedPstimeS.O(L2$bbov)
plot(p)

getMaximaLocs(L2$bbov$den$x, L2$bbov$den$y) ## take the first one
L2$bbov <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[2]], L1$bbov$cell.cycle.genes.df, 
                                           L1$bbov$sds.data, bd.tc.logCPM, lag.time = 29)


p <- plotUpdatedPstimeS.O(L2$bbov)

plot(p)

## BDIV_C
L1$bdiv_cow <- fitPseudoTime(S.O.integrated.list[[3]], reverset.time = T)
L2$bdiv_cow <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[3]], L1$bdiv_cow$cell.cycle.genes.df, 
                                               L1$bdiv_cow$sds.data, bd.tc.logCPM)

plot(L2$bdiv_cow$den)
print(L2$bdiv_cow$lag.time)
p <- plotUpdatedPstimeS.O(L2$bdiv_cow)
plot(p)

## Adjusting the start
getMaximaLocs(L2$bdiv_cow$den$x, L2$bdiv_cow$den$y) ## take the first one

L2$bdiv_cow <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[3]], L1$bdiv_cow$cell.cycle.genes.df, 
                                               L1$bdiv_cow$sds.data, bd.tc.logCPM, lag.time = 29)

p <- plotUpdatedPstimeS.O(L2$bdiv_cow)

plot(p)


## BDIV_H
L1$bdiv_human <- fitPseudoTime(S.O.integrated.list[[4]], reverset.time = T)
L2$bdiv_human <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[4]], L1$bdiv_human$cell.cycle.genes.df, 
                                                 L1$bdiv_human$sds.data, bd.tc.logCPM)

plot(L2$bdiv_human$den)
print(L2$bdiv_human$lag.time)

p <- plotUpdatedPstimeS.O(L2$bdiv_human)
plot(p)

## Adjusting the start
getMaximaLocs(L2$bdiv_human$den$x, L2$bdiv_human$den$y) ## take the first one

L2$bdiv_human <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[4]], L1$bdiv_human$cell.cycle.genes.df, 
                                                 L1$bdiv_human$sds.data, bd.tc.logCPM, lag.time = 28)

p <- plotUpdatedPstimeS.O(L2$bdiv_human)
plot(p)


## BMIC
L1$bmic <- fitPseudoTime(S.O.integrated.list[[5]], reverset.time = T)
L2$bmic <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[5]], L1$bmic$cell.cycle.genes.df, 
                                           L1$bmic$sds.data, bd.tc.logCPM)

plot(L2$bmic$den)
print(L2$bmic$lag.time)
p <- plotUpdatedPstimeS.O(L2$bmic)
plot(p)

## Adjusting the start
getMaximaLocs(L2$bmic$den$x, L2$bmic$den$y) ## take the first one


L2$bmic <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[5]], L1$bmic$cell.cycle.genes.df, 
                                           L1$bmic$sds.data, bd.tc.logCPM, lag.time = 31)

p <- plotUpdatedPstimeS.O(L2$bmic)
plot(p)


saveRDS(L1, '../Input/compScBabesia/RData/all_pstime_fits_L1.RData')
saveRDS(L2, '../Input/compScBabesia/RData/all_pstime_fits_L2.RData')


cell_cycle_genes_df <- lapply(L1, `[[`, 1)
saveRDS(cell_cycle_genes_df, '../Input/compScBabesia/RData/all_cell_cycle_genes_df.RData')

cell_cycle_genes_df_adj <- lapply(L2, `[[`, 8)
saveRDS(cell_cycle_genes_df_adj, '../Input/compScBabesia/RData/all_cell_cycle_genes_df_adj.RData')

sc.tc.df.adj <- lapply(L2, `[[`, 10)
saveRDS(sc.tc.df.adj, '../Input/compScBabesia/RData/all_sc_tc_df_adj.RData')


#### No Rerun above

titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)", "B. mic")
ps <- lapply(1:length(titles), function(i){
  p <- ggplot(L2[[i]]$pc.sds.adj, aes(x=-PC_1,y=PC_2)) +
    geom_point(aes(
      fill = cluster.y
    ), shape=21, size = 1.5)+ 
    geom_path(aes(x=sc1[cell.ord],y=sc2[cell.ord])) +
    geom_point(aes(x=-sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'black', shape=21, size = 5, stroke = 1.2)+
    geom_point(aes(x=-sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'blue', shape=8, size = 4, stroke = 1.1)+
    theme_bw(base_size = 14) +
    theme(legend.position = "none") +
    ylab('PC2') + xlab('PC1') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(titles[i]) +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) +
    guides(color = 'none')
  
})


p <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]],ps[[5]], ncol = 5)

plot(p)

ggsave(filename="../Output/compScBabsSpecies/figs/pstime_fits.png",
       plot=p,
       width = 12, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


#LabelClusters(plot = p1, id = 'adj.time.idx', size = 7, repel = T, color = 'black')


## Generate gam genes S.O update an apply network smoothing on common genes
gam.genes <- lapply(L1, function(L){
  gam.genes <- L$gam.genes
})

gam.common <- Reduce(intersect, gam.genes)

S.O.integrated.list.update.orth <- lapply(1:length(L2), function(i){
  S.O <- subset(L2[[i]]$S.O.bd.update, features = gam.common)
  
  return(subset(L2[[i]]$S.O.bd.update, features = gam.common))
})

names(S.O.integrated.list.update.orth) <- names(S.O.list)
saveRDS(S.O.integrated.list.update.orth, '../Input/compScBabesia/RData/S.O.integrated.list.pstime.GAM.intersect.Rdata')

S.O.integrated.list.update <- lapply(L2, `[[`, 1)

names(S.O.integrated.list.update) <- names(S.O.list)
saveRDS(S.O.integrated.list.update, '../Input/compScBabesia/RData/S.O.integrated.list.pstime.Rdata')


S.O.integrated.list.update.orth.indiv <- lapply(1:length(L2), function(i){
  S.O <- subset(L2[[i]]$S.O.bd.update, features = gam.genes[[i]])
  
  return(S.O)
})

names(S.O.integrated.list.update.orth.indiv) <- names(S.O.list)
saveRDS(S.O.integrated.list.update.orth.indiv, '../Input/compScBabesia/RData/S.O.integrated.list.pstime.GAM.indiv.Rdata')

# S.O.gams <- lapply(S.O.integrated.list.update, function(S.O){
#   S.O.gam <- subset(S.O, features = gam.common)
#   #S.O.gam <- prep_S.O(S.O.gam)
#   #S.O.gam.smooth <- smooth.S.O(S.O.gam)
#   return(S.O.gam)
# })

#S.O.gam = lapply(S.O.gams, `[[`, 1)
#S.O.gam.smooth = lapply(S.O.gams, `[[`, 2)

#saveRDS(S.O.gams, '../Input/compScBabesia/RData/S.O.list.ortholog.gam.Rdata')
#saveRDS(S.O.gam.smooth, '../Input/compScBabesia/RData/S.O.list.ortholog.gam.smooth.Rdata')


######
Assays(S.O.gams[[1]])


my.gene <- 'Bdiv-015780c' ## ASF

DefaultAssay(S.O.gams[[2]]) <- 'RNA'
p2 <- FeaturePlot(object = S.O.gams[[5]], 
                  #shape.by = 'spp',
                  #split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = my.gene,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)


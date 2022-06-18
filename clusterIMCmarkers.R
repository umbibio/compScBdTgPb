library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(cowplot)
library(gam)
library(princurve)
library(parallel)
library(sme)
library(plotly)
library(bigmemory)
library(dtwclust)
library(doParallel)



source('./util_funcs.R')


S.O <- readRDS('../Input/compScBdTgPb/RData/S.O.toxo_MJ_lables.RData')
sc.rna.cell.cycle.genes.df.adj <- readRDS('../Output/scClockOut/forDavid/rna_seq_gene_expr_pseudo_time.RData')
sc.atac.cell.cycle.genes.df.adj <- readRDS('../Output/scClockOut/forDavid/atac_seq_gene_access_pseudo_time.RData')
IMC_markers <- read.xlsx('../Input/compScBdTgPb/gene_function/IMC genes Toxoplasma.xlsx', sheet = 1)
ME49_TGGT1  <- read.xlsx('../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx')

IMC_markers <- left_join(IMC_markers, ME49_TGGT1, by = c('Gene.ID' = 'TGGT1'))

#comm.genes <- intersect(unique(sc.rna.cell.cycle.genes.df.adj$GeneID), unique(sc.rna.cell.cycle.genes.df.adj$GeneID))

sc.rna.genes <- IMC_markers$TGME49[which(IMC_markers$TGME49 %in% unique(sc.rna.cell.cycle.genes.df.adj$GeneID))]
sc.atac.genes <- IMC_markers$TGME49[which(IMC_markers$TGME49 %in% unique(sc.atac.cell.cycle.genes.df.adj$GeneID))]

sc.rna.spline.fits <- lapply(1:length(sc.rna.genes), function(i){
  tmp <- sc.rna.cell.cycle.genes.df.adj %>% dplyr::filter(GeneID == sc.rna.genes[i]) %>%
    transmute(GeneID = GeneID, x = adj.pt, y = log2.expr)
  sc.rna.sp <- smooth.spline(tmp$x, tmp$y, lambda = 1)
  sc.rna.sp <- predict(sc.rna.sp, seq(0, 6, by = 1/20)) 
  mu <- data.frame(x = sc.rna.sp$x, y = scale(sc.rna.sp$y))
  #mu <- data.frame(x = sc.rna.sp$x, y = sc.rna.sp$y)
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
})

sc.rna.spline.fits <- bind_rows(sc.rna.spline.fits)

sc.atac.spline.fits <- lapply(1:length(sc.atac.genes), function(i){
  tmp <- sc.atac.cell.cycle.genes.df.adj %>% dplyr::filter(GeneID == sc.atac.genes[i]) %>%
    transmute(GeneID = GeneID, x = adj.pt, y = log2.expr)
  sc.atac.sp <- smooth.spline(tmp$x, tmp$y, lambda = 1)
  sc.atac.sp <- predict(sc.atac.sp, seq(0, 6, by = 1/20)) 
  mu <- data.frame(x = sc.atac.sp$x, y = scale(sc.atac.sp$y))
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
})

sc.atac.spline.fits <- bind_rows(sc.atac.spline.fits)



## Clustering

sc.rna.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
  as.data.frame()

sc.atac.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
  as.data.frame()




# sc.rna.IMC <- sc.rna.wide[,colnames(sc.rna.wide) %in% IMC_markers$TGME49]
# sc.atac.IMC <- sc.atac.wide[,colnames(sc.atac.wide) %in% IMC_markers$TGME49]

sc.rna.IMC <- sc.rna.wide[,2:ncol(sc.rna.wide)]
sc.atac.IMC <- sc.atac.wide[,2:ncol(sc.atac.wide)]


sc.rna.IMC.markers.hc_dtw <- dtwClustCurves(sc.rna.IMC, nclust = 4L)
sc.atac.IMC.markers.hc_dtw <- dtwClustCurves(sc.atac.IMC, nclust = 4L)

plot(sc.rna.IMC.markers.hc_dtw, type = 'sc')
plot(sc.atac.IMC.markers.hc_dtw, type = 'sc')


## GGplot cluster graphs

sc.rna.long <- sc.rna.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'y') %>%
  as.data.frame()

sc.rna.long <- left_join(sc.rna.long, IMC_markers, by = c('GeneID' = 'TGME49'))
colnames(sc.rna.long) <- c('time', 'GeneID', 'normExpr', 'TGGT1', 'Name')

sc.atac.long <- sc.atac.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'y') %>%
  as.data.frame()

sc.atac.long <- left_join(sc.atac.long, IMC_markers, by = c('GeneID' = 'TGME49'))
colnames(sc.atac.long) <- c('time', 'GeneID', 'normExpr', 'TGGT1', 'Name')

sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna.IMC), cluster = cutree(sc.rna.IMC.markers.hc_dtw, k = 4))
sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac.IMC), cluster = cutree(sc.atac.IMC.markers.hc_dtw, k = 4))

sc.rna.long.clust <- inner_join(sc.rna.long, sc.rna.clust.info, by = 'GeneID')
sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')


p1 <- lapply(unique(sc.rna.long.clust$cluster), function(i){
  tmp <- sc.rna.long.clust %>% dplyr::filter(cluster == i)
  ggplot(tmp, aes(x = time, y = normExpr, color = Name)) + 
    geom_line() + theme_bw()
})


p <- do.call(grid.arrange, c(p1, nrow=1))

ggsave(p, file = "../Output/compScBdTgPb/figs/toxo_IMC_markers_cluster_smoothed_v2.pdf.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 6) # The height of the plot in inches





p2 <- lapply(unique(sc.atac.long.clust$cluster), function(i){
  tmp <- sc.atac.long.clust %>% dplyr::filter(cluster == i)
  ggplot(tmp, aes(x = time, y = normExpr, color = Name)) + 
    geom_line() + theme_bw()
})


p <- do.call(grid.arrange, c(p2, nrow=2))

ggsave(p, file = "../Output/compScBdTgPb/figs/toxo_atac_IMC_markers_cluster_smoothed.pdf",   # The directory you want to save the file in
       width = 8, # The width of the plot in inches
       height = 8) # The height of the plot in inches

IMC_markers

FeaturePlot(S.O, features = 'TGGT1-253470', reduction = 'pca', label = T)

pdf(file = "../Output/scClockFigs/sc_egress_markers_expression_curve_clusters.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches





##
sc.egress.markers.tc.mus.long <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df.adj$variable), extend = F)
colnames(sc.egress.markers.tc.mus.long) <- c('GeneID', 't', 'y')
sc.egress.markers.tc.mus.long <- sc.egress.markers.tc.mus.long[sc.egress.markers.tc.mus.long$GeneID %in% egress_markers$GeneID,]
#sc.egress.markers.tc.mus.long <- sc.egress.markers.tc.mus.long[sc.egress.markers.tc.mus.long$GeneID %in% egress_markers.putative$BdGene,]
sc.egress.markers.tc.mus.long.stats <- sc.egress.markers.tc.mus.long %>% group_by(GeneID) %>% summarise(mean = mean(y), sd = sd(y))

sc.egress.markers.tc.mus.long <- left_join(sc.egress.markers.tc.mus.long, sc.egress.markers.tc.mus.long.stats, by = 'GeneID') %>% 
  rowwise() %>% mutate(y.st = (y - mean)/sd)
clust.info <- data.frame(GeneID = colnames(sc.egress.markers.tc.mus), cluster = cutree(sc.egress.markers.hc_dtw, k = 5))

sc.egress.markers.tc.mus.long <- inner_join(sc.egress.markers.tc.mus.long, clust.info, by = 'GeneID')
sc.egress.markers.tc.mus.long$Name <- egress_markers$Bd_name[match(sc.egress.markers.tc.mus.long$GeneID, egress_markers$GeneID)]


p <- lapply(unique(sc.egress.markers.tc.mus.long$cluster), function(i){
  tmp <- sc.egress.markers.tc.mus.long %>% dplyr::filter(cluster == i)
  ggplot(tmp, aes(x = t, y = y, color = Name)) + 
    geom_line() + theme_bw()
})

do.call(grid.arrange, c(p, nrow=2))


pdf(file = "../Output/scClockFigs/sc_egress_markers_expression_curve_clusters.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches




## Generate the clusters
num.clust <- 4L


## Generate the clusters
num.clust <- 4L

sc.rna.hc_dtws <- dtwClustCurves(sc.rna.dtw.wide[2:ncol(sc.rna.dtw.wide)], nclust = num.clust)
plot(sc.rna.hc_dtws, type = 'sc')
plot(sc.rna.hc_dtws, type = "series", clus = 1L)
plot(sc.rna.hc_dtws, type = "centroids", clus = 1L)


sc.atac.hc_dtws <- dtwClustCurves(sc.atac.dtw.wide[2:ncol(sc.atac.dtw.wide)], nclust = num.clust)
plot(sc.atac.hc_dtws, type = 'sc')
plot(sc.atac.hc_dtws, type = "series", clus = 1L)
plot(sc.atac.hc_dtws, type = "centroids", clus = 1L)

sc.hc_dtw <- dtwClustCurves(sc.tc.mus[,2:ncol(sc.tc.mus)], nclust = num.clust)


sc.tc.mus.scale <- sc.tc.mus
sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                   center = T,scale = T)
sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(!starts_with('t'), names_to = 'GeneID', values_to = 'y')



sc.hc_dtw.df <- data.frame(GeneID = unique(sc.tc.mus.scale$GeneID), 
                           order = as.numeric(sc.hc_dtw$order),
                           cluster = cutree(sc.hc_dtw,k = num.clust))

sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.hc_dtw.df, by = 'GeneID')

sc.tc.mus.scale$GeneID <- factor(as.character(sc.tc.mus.scale$GeneID),
                                 levels = unique(sc.tc.mus.scale$GeneID[sc.tc.mus.scale$order]))



hc_eucledian.df <- withinCalssReOrder(sc.tc.mus.scale) 
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, hc_eucledian.df, by = 'GeneID')
sc.tc.mus.scale <- sc.tc.mus.scale %>% mutate(cluster = as.factor(cluster),
                                              GeneID.reord = reorder_within(GeneID, hc_eucledian.order, cluster)) 

## map the clusteres
sc.overlap <- matchClustersToPhase(sc.hc_dtw.df, BD.markers.sig)

sc.phase.match <- sc.overlap %>% group_by(cluster) %>% summarize(phase = markers[which.max(percent)])
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.phase.match, by = 'cluster')

## mark location in heatmap
p <- ggplot(sc.tc.mus.scale, aes(x = t, y = GeneID, fill = y)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("time/cells") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    #axis.ticks = element_blank(),
    #axis.text.y  = element_blank(),
    legend.position = "none") +
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  scale_y_discrete(breaks = egress_markers$GeneID, labels = egress_markers$Bd_name, guide=guide_axis(n.dodge=3)) + 
  #scale_y_discrete(breaks = egress_markers.putative$BdGene, labels = egress_markers.putative$BdGene, guide=guide_axis(n.dodge=3)) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=5, face="bold", color = 'black')
  ) 

plot(p)  


ggsave(filename="../Output/scClockFigs/curve_cluster_heatmap_sc_egress.pdf",
       plot=p6,
       width = 6, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)




## Recluster the egress markers

sc.egress.markers.tc.mus <- sc.tc.mus[,colnames(sc.tc.mus) %in% egress_markers$GeneID]
#sc.egress.markers.tc.mus <- sc.tc.mus[,colnames(sc.tc.mus) %in% egress_markers.putative$BdGene]

sc.egress.markers.hc_dtw <- dtwClustCurves(sc.egress.markers.tc.mus, nclust = 6L)

plot(sc.egress.markers.hc_dtw, type = 'sc')


##
sc.egress.markers.tc.mus.long <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df.adj$variable), extend = F)
colnames(sc.egress.markers.tc.mus.long) <- c('GeneID', 't', 'y')
sc.egress.markers.tc.mus.long <- sc.egress.markers.tc.mus.long[sc.egress.markers.tc.mus.long$GeneID %in% egress_markers$GeneID,]
#sc.egress.markers.tc.mus.long <- sc.egress.markers.tc.mus.long[sc.egress.markers.tc.mus.long$GeneID %in% egress_markers.putative$BdGene,]
sc.egress.markers.tc.mus.long.stats <- sc.egress.markers.tc.mus.long %>% group_by(GeneID) %>% summarise(mean = mean(y), sd = sd(y))

sc.egress.markers.tc.mus.long <- left_join(sc.egress.markers.tc.mus.long, sc.egress.markers.tc.mus.long.stats, by = 'GeneID') %>% 
  rowwise() %>% mutate(y.st = (y - mean)/sd)
clust.info <- data.frame(GeneID = colnames(sc.egress.markers.tc.mus), cluster = cutree(sc.egress.markers.hc_dtw, k = 5))

sc.egress.markers.tc.mus.long <- inner_join(sc.egress.markers.tc.mus.long, clust.info, by = 'GeneID')
sc.egress.markers.tc.mus.long$Name <- egress_markers$Bd_name[match(sc.egress.markers.tc.mus.long$GeneID, egress_markers$GeneID)]


p <- lapply(unique(sc.egress.markers.tc.mus.long$cluster), function(i){
  tmp <- sc.egress.markers.tc.mus.long %>% dplyr::filter(cluster == i)
  ggplot(tmp, aes(x = t, y = y, color = Name)) + 
    geom_line() + theme_bw()
})

do.call(grid.arrange, c(p, nrow=2))


pdf(file = "../Output/scClockFigs/sc_egress_markers_expression_curve_clusters.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches

par(mfrow = c(1,2))
v <- CDPK4.id
ind <- which(unique(sc.tc.df.adj$variable) == v)
plot.sme(sc.tc.fits[[ind]], paste('sc', v))
ind <- which(unique(sync.tc.df$variable) == v)
plot.sme(sync.tc.fits[[ind]], paste('sync', v), conf = F)

dev.off()







#################

prod.desc[grep('calcium',prod.desc$Product.Description),]
CDPK4.id <- "Bdiv_024410"
prod.desc[grep('GMP',prod.desc$Product.Description),]
PKG.id <- "Bdiv_020500"

sc.tc.mus <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df.adj$variable), extend = F)
colnames(sc.tc.mus) <- c('GeneID', 't', 'y')
sc.tc.mus <- sc.tc.mus %>% dplyr::filter(GeneID %in% c(BD.markers.sig$GeneID, PKG.id) )%>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  as.data.frame()

## Generate the clusters
num.clust <- 4L


sc.hc_dtw <- dtwClustCurves(sc.tc.mus, nclust = num.clust)


sc.tc.mus.scale <- sc.tc.mus
sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                   center = T,scale = T)
sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(!starts_with('t'), names_to = 'GeneID', values_to = 'y')



sc.hc_dtw.df <- data.frame(GeneID = unique(sc.tc.mus.scale$GeneID), 
                           order = as.numeric(sc.hc_dtw$order),
                           cluster = cutree(sc.hc_dtw,k = num.clust))

sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.hc_dtw.df, by = 'GeneID')

sc.tc.mus.scale$GeneID <- factor(as.character(sc.tc.mus.scale$GeneID),
                                 levels = unique(sc.tc.mus.scale$GeneID[sc.tc.mus.scale$order]))



hc_eucledian.df <- withinCalssReOrder(sc.tc.mus.scale) 
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, hc_eucledian.df, by = 'GeneID')
sc.tc.mus.scale <- sc.tc.mus.scale %>% mutate(cluster = as.factor(cluster),
                                              GeneID.reord = reorder_within(GeneID, hc_eucledian.order, cluster)) 

## map the clusteres
sc.overlap <- matchClustersToPhase(sc.hc_dtw.df, BD.markers.sig)

sc.phase.match <- sc.overlap %>% group_by(cluster) %>% summarize(phase = markers[which.max(percent)])
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.phase.match, by = 'cluster')



p4 <- ggplot(sc.tc.mus.scale, aes(x = t, y = GeneID, fill = y)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("time/cells") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    #axis.ticks = element_blank(),
    #axis.text.y  = element_blank(),
    legend.position = "none") +
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=8, face="bold")
  ) + scale_y_discrete(breaks = c(CDPK4.id, PKG.id), labels = c('CDPK4', 'PKG'))

plot(p4)  


ggsave(filename="../Output/scClockFigs/curve_cluster_heatmap_sc_CDPK4_PKG.pdf",
       plot=p4,
       width = 5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)

## CDPK4 trend

pdf(file = "../Output/scClockFigs/CDPK4_expression_curve.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches

par(mfrow = c(1,2))
v <- CDPK4.id
ind <- which(unique(sc.tc.df.adj$variable) == v)
plot.sme(sc.tc.fits[[ind]], paste('sc', v))
ind <- which(unique(sync.tc.df$variable) == v)
plot.sme(sync.tc.fits[[ind]], paste('sync', v), conf = F)

dev.off()

pdf(file = "../Output/scClockFigs/PKG_expression_curve.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches


par(mfrow = c(1,2))
v <- PKG.id
ind <- which(unique(sc.tc.df.adj$variable) == v)
plot.sme(sc.tc.fits[[ind]], paste('sc', v))
ind <- which(unique(sync.tc.df$variable) == v)
plot.sme(sync.tc.fits[[ind]], paste('sync', v), conf = F)

dev.off()


#### Expression
Idents(S.O.bd.filt) <- 'seurat_clusters'
DefaultAssay(S.O.bd.filt) <- "RNA"


p2 <- FeaturePlot(object = S.O.bd.filt, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)

p3 <- FeaturePlot(object = S.O.bd.filt, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', PKG.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p3)



ggsave(filename="../Output/scClockFigs/expression_CDPK4.pdf",
       plot=p2,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

ggsave(filename="../Output/scClockFigs/expression_PKG.pdf",
       plot=p3,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



####

p6 <- ggplot(sc.tc.mus.scale, aes(x = t, y = GeneID, fill = y)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("time/cells") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    #axis.ticks = element_blank(),
    #axis.text.y  = element_blank(),
    legend.position = "none") +
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  scale_y_discrete(breaks = egress_markers$GeneID, labels = egress_markers$Pf_NAMES, guide=guide_axis(n.dodge=3)) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=5, face="bold", color = 'black')
  ) 

plot(p6)  


ggsave(filename="../Output/scClockFigs/curve_cluster_heatmap_sc_egress.pdf",
       plot=p6,
       width = 6, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)


###
S.O.bd.filt@meta.data$Sample <- rownames(S.O.bd.filt@meta.data)
PKG.expressing.ind <- which(S.O.bd.filt[["RNA"]]@data[grep(gsub('_', '-', PKG.id), rownames(S.O.bd.filt[["RNA"]]@data)),] > log2(3))
PKG.expressing <- colnames(S.O.bd.filt[["RNA"]]@data)[PKG.expressing.ind]
S.O.bd.filt@meta.data$PKG <- ifelse(S.O.bd.filt@meta.data$Sample %in% PKG.expressing, 1, 0)
S.O.bd.filt@meta.data$PKG_clust <- paste(S.O.bd.filt@meta.data$PKG, S.O.bd.filt@meta.data$seurat_clusters, sep = '_') 
Idents(S.O.bd.filt) <- 'PKG_clust'
PKG.markers <- FindAllMarkers(object = S.O.bd.filt, only.pos = TRUE, min.pct = 0)

PKG.markers$GeneID <- gsub('-', '_', PKG.markers$gene)
PKG.markers.top <- PKG.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
PKG.markers.sig <- PKG.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05 & cluster == 1)
PKG.markers.sig <- left_join(PKG.markers.sig, prod.desc, by = 'GeneID')
write.xlsx(PKG.markers.sig, '../Output/scClockOut/PKG_expressing_cells_markers.xlsx')
FeaturePlot(object = S.O.bd.filt, 
            features = PKG.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")


WhichCells(S.O.bd.filt, slot = 'data', expression = gsub('_', '-', PKG.id) > 0 )

subset(S.O.bd.filt, subset = gsub('_', '-', PKG.id) > 1)

xx<-colnames(S.O.bd.filt)
pc.sds.adj <- readRDS('../Input/scClock/pc.sds.adj.RData')

pc.sds.adj$seurat_clusters <- S.O.bd.filt@meta.data$seurat_clusters[match(pc.sds.adj$Sample, S.O.bd.filt@meta.data$Sample)]
pc.sds.adj <- pc.sds.adj %>% mutate(phase = ifelse(seurat_clusters == 1, 'S/M', 
                                                   ifelse(seurat_clusters == 0, 'G1a',
                                                          ifelse(seurat_clusters == 3, 'G1b','G1c'))))
p <- ggplot(pc.sds.adj, aes(x=PC_1,y=PC_2)) +
  geom_point(aes(
    fill = phase
  ), shape=21, size = 1.5)+
  geom_path(aes(x=sc1[cell.ord],y=sc2[cell.ord])) +
  geom_point(aes(x=sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'black', shape=21, size = 5, stroke = 1.2)+
  geom_point(aes(x=sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'blue', shape=8, size = 4, stroke = 1.1)+
  theme_bw(base_size = 14) +
  # theme(legend.direction="horizontal",
  #       legend.position = c(0.6,0.99)) +
  theme(legend.position = "right") +
  
  ylab('PC2') + xlab('PC1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  guides(color = FALSE)

plot(p)


ggsave(filename="../Output/scClockFigs/PCA_time.pdf",
       plot=p,
       width = 6, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

##### Matched marker analysis


S.O.bd.filt@meta.data$Sample <- rownames(S.O.bd.filt@meta.data)
PKG.expressing.ind <- which(S.O.bd.filt[["RNA"]]@data[grep(gsub('_', '-', PKG.id), rownames(S.O.bd.filt[["RNA"]]@data)),] > 0)
PKG.expressing <- colnames(S.O.bd.filt[["RNA"]]@data)[PKG.expressing.ind]
S.O.bd.filt@meta.data$PKG <- ifelse(S.O.bd.filt@meta.data$Sample %in% PKG.expressing, 1, 0)
S.O.bd.filt@meta.data$PKG_clust <- paste(S.O.bd.filt@meta.data$PKG, S.O.bd.filt@meta.data$seurat_clusters, sep = '_') 
Idents(S.O.bd.filt) <- 'PKG_clust'


objs <- unique(S.O.bd.filt@meta.data$PKG)
ref.obj <- objs[1]
query.obj <- objs[2]
clusters <- unique(S.O.bd.filt@meta.data$seurat_clusters)
ref.clusts <- data.frame(ref = paste(ref.obj, clusters, sep = '_'))
ref.clusts$cluster <- gsub('.*_', '', ref.clusts$ref)
ref.clusts$dummy <- 1
query.clusts <- data.frame(query = paste(rep(query.obj, each = length(clusters)), clusters, sep = '_'))
query.clusts$dummy <- 1
query.clusts$cluster <- gsub('.*_', '', query.clusts$query)

contrasts <- full_join(ref.clusts, query.clusts, by = 'dummy') %>% 
  dplyr::filter(cluster.x == cluster.y) %>%
  transmute(ref = ref, query = query)

contrasts$cluster <- gsub('.*_', '', contrasts$ref)

#ident.1 case, ident.2 is control
Idents(S.O.bd.filt) <- 'PKG_clust'
matched.DEGs <- mclapply(split(contrasts, seq(nrow(contrasts))), function(x){
  tmp <- FindMarkers(S.O.bd.filt, ident.1 = x$query, ident.2 = x$ref, verbose = T)
  tmp$ref <- x$ref
  tmp$query <- x$query
  tmp$cluster <- x$cluster
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})


matched.DEGs <- bind_rows(matched.DEGs)
PKG.markers <- matched.DEGs %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
PKG.markers$GeneID <- gsub('-', '_', PKG.markers$gene)
PKG.markers <- left_join(PKG.markers, prod.desc, by = 'GeneID')
write.xlsx(PKG.markers, '../Output/scClockOut/PKG_expressing_cells_markers_matched_clusters.xlsx')
FeaturePlot(object = S.O.bd.filt, 
            features = unique(PKG.markers$gene), 
            cols = c("grey", "blue"), reduction = "pca")

###################3
## Do a network smoothing
## Smoothing over KNN
S.O.bd.filt <- FindNeighbors(S.O.bd.filt, reduction = "pca", dims = 1:10)

AdjacencyMat <- as.matrix(S.O.bd.filt@graphs$RNA_nn)

con <- colSums(AdjacencyMat)
con <- unlist(lapply(con, function(x) ifelse(x == 0, 0, 1/x)))

## Scaling adjacancy
for(i in 1:ncol(AdjacencyMat)){
  AdjacencyMat[,i] <- con[i] * AdjacencyMat[,i]
}

## Get the expression data
expr.norm <- as.matrix(S.O.bd.filt[["RNA"]]@data)
expr.norm.filt <- expr.norm[which(rownames(expr.norm) %in% gsub('_', '-', c(CDPK4.id, PKG.id))), ]

## Smoothing for a few iterations
max.smoothing <- 20
alpha <- 0.7
scDat <- expr.norm.filt
Ft <- scDat
for(i in 1:max.smoothing){
  Ft <- alpha * Ft %*% AdjacencyMat + (1 - alpha) * scDat
}

S.O.smooth <- S.O.bd.filt
smooth_assay <- CreateAssayObject(counts = Ft)
S.O.smooth[["smooth"]] <- smooth_assay


Assays(S.O.smooth)
Idents(S.O.smooth) <- 'seurat_clusters'

DefaultAssay(S.O.smooth) <- "RNA"

p3 <- FeaturePlot(object = S.O.smooth, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p3)



DefaultAssay(S.O.smooth) <- "smooth"

p4 <- FeaturePlot(object = S.O.smooth, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p4)





#####
S.O.ave <- S.O.bd.filt
ave_assay <- CreateAssayObject(counts = scDat)
S.O.ave[["ave"]] <- ave_assay


Assays(S.O.ave)
Idents(S.O.ave) <- 'seurat_clusters'

DefaultAssay(S.O.ave) <- "RNA"

p3 <- FeaturePlot(object = S.O.ave, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p3)

DefaultAssay(S.O.ave) <- "ave"

p4 <- FeaturePlot(object = S.O.ave, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p4)


library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(slingshot)
library(gam)
library(princurve)
library(parallel)
library(dtwclust)
library(doParallel)
library(tidyverse)
library(tidytext)
library(fda)
library(sme)


#library(sctransform)

source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
## Fit a pseudo-time curve and align using sync data
S.O.tg <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.RData')
pc.tg <- getPCA(S.O.tg)
sds.data <- getPrinCurve(pc.tg)
pc.sds.tg <- left_join(pc.tg, sds.data, by = "Sample")
pc.sds.tg$phase <- S.O.tg@meta.data$phase[match(pc.sds.tg$Sample, S.O.tg@meta.data$Sample)]
pc.sds.tg$phase <- factor(pc.sds.tg$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))

saveRDS(pc.sds.tg, '../Input/compScBdTgPb/RData/pc.sds.tg.RData')

p <- ggplot(pc.sds.tg, aes(x=PC_1,y=PC_2)) + 
  geom_point(aes(
    fill = phase
  ), shape=21, size = 1.5)+
  geom_path(aes(x=sc1[cell.ord],y=sc2[cell.ord])) + 
  theme_bw(base_size = 14) + 
  theme(legend.position=c(1,1),legend.justification=c(1,1), 
        legend.title = element_blank(),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.direction="vertical") + 
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

ggsave(filename="../Output/compScBdTgPb/figs/initial_pseudo_time_TG2.pdf",
       plot=p,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)




## Pseudo-time analysis with SLingshot
## Identify genes that correlate with time

#Y <- log2(S.O.bd.filt@assays$smooth@counts + 1) ## smoothed version
Y <- log2(S.O.tg@assays$RNA@counts + 1)
var.genes <- names(sort(apply(Y, 1, var),decreasing = TRUE))#[1:1000] 
Y <- Y[var.genes, ]

pt <- sds.data$pt

## Map the pseudo-time to 0-12:20 hours 
#t <- (12 + 1/3) * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
t <- (6 + 1/6) * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))

sds.data$t <- t

## time-index cells in 20 min intervals and identify cells in each partition
## They will be considered as replicates
#time.breaks <- seq(1/3, 12 + 1/3, by = 1/3) 
time.breaks <- seq(1/6, 6 + 1/6, by = 1/6) 
time.idx <- rep(0, nrow(sds.data))

ind <- which(sds.data$t <= time.breaks[1])
time.idx[ind] <- 0

for(i in 2:length(time.breaks)){
  ind <- which(sds.data$t > time.breaks[(i-1)] & sds.data$t <= time.breaks[i])
  time.idx[ind] <- i - 1
}

sds.data$time.idx <- time.idx

## Update the time to 20 min increments
#sds.data$t <- (time.idx) * (1/3)
sds.data$t <- (time.idx) * (1/6)

sds.data <- sds.data %>%  
  group_by(time.idx) %>% mutate(rep = seq(1:n()))


rownames(sds.data) <- sds.data$Sample

## Add the new clusters as meta-data
S.O.tg <- AddMetaData(S.O.tg, sds.data)
saveRDS(S.O.tg, '../Input/compScBdTgPb/RData/S.O.tg.RData')


## Run a GAM regression of expression on the pseudo-time
## Use parallel computation to speed things up. 16 cores
gam.pval <- mclapply(1:nrow(Y), function(z){
  d <- data.frame(z=as.numeric(Y[z,]), t=as.numeric(pt))
  tmp <- gam(z ~ lo(t), data=d)
  #p <- summary(tmp)[4][[1]][1,5]
  p <- summary(tmp)$anova$`Pr(F)`[2]
  p
}, mc.cores = num.cores)

gam.pval <- unlist(gam.pval)
names(gam.pval) <- rownames(Y)
## Remove the NA's and get the best fits
#gam.pval <- gam.pval[-which(is.na(gam.pval))]
gam.pval.adj <- p.adjust(gam.pval, method = 'fdr', n = length(gam.pval))
gam.pval.sig <- gam.pval[gam.pval.adj < 0.01] 
print(length(gam.pval.sig)) ## number of correlating genes

## Sort the cells on the pt
cell.ord <- sds.data$cell.ord

topgenes <- names(sort(gam.pval.sig, decreasing = FALSE))  
cell.cycle.genes.expr <- as.matrix(S.O.tg@assays$RNA@data[topgenes, cell.ord])
#cell.cycle.genes.expr <- as.matrix(S.O.bd.filt@assays$smooth@data[topgenes, cell.ord]) ## smoothed version


cell.cycle.genes.df <- data.frame(GeneID = rownames(cell.cycle.genes.expr),
                                  cell.cycle.genes.expr) %>% 
  pivot_longer(-c(GeneID), names_to = 'Sample', values_to = 'log2.expr')


cell.cycle.genes.df$GeneID <- gsub('-', '_', cell.cycle.genes.df$GeneID)
cell.cycle.genes.df <- left_join(cell.cycle.genes.df, sds.data, by = 'Sample')
cell.cycle.genes.df$cluster <- S.O.tg@meta.data$seurat_clusters[match(cell.cycle.genes.df$Sample, 
                                                                           rownames(S.O.tg@meta.data))]

cell.cycle.genes.df$phase <- S.O.tg@meta.data$phase[match(cell.cycle.genes.df$Sample, 
                                                                           rownames(S.O.tg@meta.data))]


saveRDS(cell.cycle.genes.df, '../Input/compScBdTgPb/RData/tg_cell_cycle_genes_df.RData')

## Filtering to include genes that fit well with pseudo time
S.O.tg.gam <- subset(S.O.tg, features = names(gam.pval.sig))
S.O.tg.gam <- prep_S.O(S.O.tg.gam)
S.O.tg.gam <- smooth.S.O(S.O.tg.gam)
saveRDS(S.O.tg.gam, '../Input/compScBdTgPb/RData/S.O.tg.gam.RData')


sc.tc.df <- cell.cycle.genes.df %>% 
  transmute(y = log2.expr, tme = t, ind = rep, variable = GeneID)
saveRDS(sc.tc.df, '../Input/compScBdTgPb/RData/tg_sc_tc_df.RData')


sc.tc.fits <- mclapply(unique(sc.tc.df$variable),
                       function(v)
                         sme(sc.tc.df[sc.tc.df$variable==v,c("y","tme","ind")],
                             lambda.mu = 6, lambda.v = 6), mc.cores = num.cores)

#saveRDS(object = sc.tc.fits, file = "../Input/compScBdTgPb/RData/tg_sme_fits_sc_tc_20min.RData")
saveRDS(object = sc.tc.fits, file = "../Input/compScBdTgPb/RData/tg_sme_fits_sc_tc_10min.RData")


## Plot a few curves to check the alignments

vs = unique(sc.tc.df$variable)[1:16]


pdf(file = "../Output/compScBdTgPb/figs/tg_sme_fits_sc.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

par(mfrow = c(4,4))
for(v in vs){
  ind <- which(unique(sc.tc.df$variable) == v)
  plot.sme(sc.tc.fits[[ind]], v)
}

dev.off()



## Get the sme mean splines and filter to include marker genes only
TG.markers.sig <- readRDS('../Input/compScBdTgPb/RData/TG.markers.sig.RData')

sc.tc.mus <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df$variable), extend = F)
colnames(sc.tc.mus) <- c('GeneID', 't', 'y')
sc.tc.mus <- sc.tc.mus %>% dplyr::filter(GeneID %in% TG.markers.sig$GeneID) %>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  as.data.frame()

## Generate the clusters
num.clust <- 8L

sc.hc_dtw <- dtwClustCurves(sc.tc.mus[,2:ncol(sc.tc.mus)], nclust = num.clust)
plot(sc.hc_dtw, type = 'sc')
plot(sc.hc_dtw, type = 'centroids')
plot(sc.hc_dtw, type = "series", clus = 2L)
plot(sc.hc_dtw, type = "centroids", clus = 2L)

saveRDS(sc.hc_dtw, '../Input/compScBdTgPb/RData/tg_sc.hc_dtw.RData')

## Scale the data for heatmaps
sc.tc.mus.scale <- sc.tc.mus
sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                   center = T,scale = T)
sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(starts_with('TG'), names_to = 'GeneID', values_to = 'y')


## Add curve cluster info

sc.hc_dtw.df <- data.frame(GeneID = unique(sc.tc.mus.scale$GeneID), 
                           order = as.numeric(sc.hc_dtw$order),
                           cluster = cutree(sc.hc_dtw,k = num.clust))

sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.hc_dtw.df, by = 'GeneID')

#sc.tc.mus.scale$GeneID <- factor(as.character(sc.tc.mus.scale$GeneID),
#                                 levels = unique(sc.tc.mus.scale$GeneID[sc.tc.mus.scale$order]))



## Reorder the genes within each cluster.
#hc_eucledian.df <- withinCalssReOrder(sc.tc.mus.scale) 
#sc.tc.mus.scale <- left_join(sc.tc.mus.scale, hc_eucledian.df, by = 'GeneID')
#sc.tc.mus.scale <- sc.tc.mus.scale %>% mutate(cluster = as.factor(cluster),
#                                              GeneID.reord = reorder_within(GeneID, hc_eucledian.order, cluster)) 

## map the clusters
sc.peak.order <- sc.tc.mus.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(t, y))
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.peak.order, by = 'GeneID')

TG.markers.sig$cluster <- gsub('G1.b', 'G1', gsub('G1.a', 'G1', as.character(TG.markers.sig$cluster)))
sc.overlap <- matchClustersToPhase(sc.hc_dtw.df, TG.markers.sig)
sc.overlap$cluster <- as.numeric(sc.overlap$cluster)

sc.phase.match <- sc.overlap %>% group_by(cluster) %>% summarize(phase = markers[which.max(percent)])
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.phase.match, by = 'cluster')

sc.tc.mus.scale$phase <- factor(sc.tc.mus.scale$phase, levels = c('G1', 'S', 'M', 'C'))

#sc.tc.mus.scale$phase <- factor(sc.tc.mus.scale$phase, levels = c('G1.a','G1.b', 'S', 'M', 'C'))


p2 <- ggplot(sc.tc.mus.scale, aes(x = t, y = reorder_within(GeneID, -peak.ord, phase), fill = y)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("time/cells") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    legend.position = "none") +
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p2)  


ggsave(filename="../Output/compScBdTgPb/figs/tg_curve_cluster_heatmap_scg_g1a_g1b.png",
       plot=p2,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)

sc.overlap$cluster <- factor(sc.overlap$cluster, levels = c(5, 4, 3, 2, 1))
sc.overlap$markers <- factor(sc.overlap$markers, levels = c('G1', 'S', 'M', 'C'))

p3 <- ggplot(sc.overlap, aes(x = cluster, y = markers, fill = percent)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("markers") + xlab("clusters") + theme_bw() + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  # theme(
  #   axis.text.x = element_blank(),
  #   axis.ticks = element_blank(),
  #   axis.text.y  = element_blank(),
  #   legend.position = "none") +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p3)  


saveRDS(sc.tc.mus, '../Input/compScBdTgPb/RData/tg_sc.tc.mus.RData')
saveRDS(sc.tc.mus.scale, '../Input/compScBdTgPb/RData/tg_sc.tc.mus.scale.RData')


### Calculate transition times

tmp <- sc.tc.mus.scale %>% ungroup() %>% group_by(phase) %>% summarise(mean.peak = median(peak.ord))

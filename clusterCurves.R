library(dtwclust)
library(doParallel)
library(tidyverse)
library(tidytext)

source('./util_funcs.R')

getCurvePeakLoc <- function(t, y){
  
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
    maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = 0.9))
    
    ## Peaks for entities of interest
    entity.x = maxval$x[maxima.outliers]
    entity.y = maxval$y[maxima.outliers]
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}

BD.markers.sig <- readRDS('../Input/compScBdTgPb/RData/BD.markers.sig.RData')

## Get the sme mean splines and filter to include marker genes only
sync.tc.df <- readRDS('../Input/compScBdTgPb/RData/bd_sync_tc_df.RData')
sync.tc.fits <- readRDS('../Input/compScBdTgPb/RData/bd_sme_fits_sync_tc_20min.RData')

sync.tc.mus <- smoothSplineSmeFits(sync.tc.fits, unique(sync.tc.df$variable), extend = T)

colnames(sync.tc.mus) <- c('GeneID', 't', 'y')
sync.tc.mus <- sync.tc.mus %>% dplyr::filter(GeneID %in% BD.markers.sig$GeneID) %>%
  pivot_wider(names_from = 'GeneID', values_from = 'y')  %>%
  as.data.frame()

sc.tc.df.adj <- readRDS('../Input/compScBdTgPb/RData/bd_sc_tc_df_adj.RData')
sc.tc.fits <- readRDS('../Input/compScBdTgPb/RData/bd_sme_fits_sc_tc_20min.RData')

sc.tc.mus <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df.adj$variable), extend = F)
colnames(sc.tc.mus) <- c('GeneID', 't', 'y')
sc.tc.mus <- sc.tc.mus %>% dplyr::filter(GeneID %in% BD.markers.sig$GeneID) %>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
  as.data.frame()

## Generate the clusters
num.clust <- 4L

sync.hc_dtw <- dtwClustCurves(sync.tc.mus[,2:ncol(sync.tc.mus)], nclust = num.clust)
plot(sync.hc_dtw, type = 'sc')
plot(sync.hc_dtw, type = "series", clus = 1L)
plot(sync.hc_dtw, type = "centroids", clus = 1L)

saveRDS(sync.hc_dtw, '../Input/compScBdTgPb/RData/bd_sync_hc_dtw.RData')

sc.hc_dtw <- dtwClustCurves(sc.tc.mus[,2:ncol(sc.tc.mus)], nclust = num.clust)
plot(sc.hc_dtw, type = 'sc')
plot(sc.hc_dtw, type = 'centroids')
plot(sc.hc_dtw, type = "series", clus = 2L)
plot(sc.hc_dtw, type = "centroids", clus = 2L)

saveRDS(sc.hc_dtw, '../Input/compScBdTgPb/RData/bd_sc_hc_dtw.RData')

## Scale the data for heatmaps
sync.tc.mus.scale <- sync.tc.mus
sync.tc.mus.scale[,2:ncol(sync.tc.mus.scale)] <- scale(sync.tc.mus.scale[,2:ncol(sync.tc.mus.scale)],
                                                       center = T,scale = T)
sync.tc.mus.scale <- sync.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(cols = -t, names_to = 'GeneID', values_to = 'y')


sc.tc.mus.scale <- sc.tc.mus
sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                       center = T,scale = T)
sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(cols = -t, names_to = 'GeneID', values_to = 'y')


## Add curve cluster info
sync.hc_dtw.df <- data.frame(GeneID = unique(sync.tc.mus.scale$GeneID), 
                        order = as.numeric(sync.hc_dtw$order),
                        cluster = cutree(sync.hc_dtw,k = num.clust))
sync.tc.mus.scale <- left_join(sync.tc.mus.scale, sync.hc_dtw.df, by = 'GeneID')

#sync.tc.mus.scale$GeneID <- factor(as.character(sync.tc.mus.scale$GeneID),
#                                   levels = unique(sync.tc.mus.scale$GeneID[sync.tc.mus.scale$order]))


sc.hc_dtw.df <- data.frame(GeneID = unique(sc.tc.mus.scale$GeneID), 
                             order = as.numeric(sc.hc_dtw$order),
                             cluster = cutree(sc.hc_dtw,k = num.clust))

sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.hc_dtw.df, by = 'GeneID')

#sc.tc.mus.scale$GeneID <- factor(as.character(sc.tc.mus.scale$GeneID),
#                                   levels = unique(sc.tc.mus.scale$GeneID[sc.tc.mus.scale$order]))



## Reorder the genes within each cluster.

sync.peak.order <- sync.tc.mus.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(t, y))
sync.tc.mus.scale <- left_join(sync.tc.mus.scale, sync.peak.order, by = 'GeneID')

sc.peak.order <- sc.tc.mus.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(t, y))
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.peak.order, by = 'GeneID')

#hc_eucledian.df <- withinCalssReOrder(sync.tc.mus.scale) 
#sync.tc.mus.scale <- left_join(sync.tc.mus.scale, hc_eucledian.df, by = 'GeneID')
#sync.tc.mus.scale <- sync.tc.mus.scale %>% mutate(cluster = as.factor(cluster),
#                                                  GeneID.reord = reorder_within(GeneID, hc_eucledian.order, cluster)) 


sc.peak.order <- sc.tc.mus.scale %>% group_by(GeneID) %>% summarise(peak.ord = which.max(y))
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.peak.order, by = 'GeneID')

#hc_eucledian.df <- withinCalssReOrder(sc.tc.mus.scale) 
#sc.tc.mus.scale <- left_join(sc.tc.mus.scale, hc_eucledian.df, by = 'GeneID')
#sc.tc.mus.scale <- sc.tc.mus.scale %>% mutate(cluster = as.factor(cluster),
#                                                 GeneID.reord = reorder_within(GeneID, hc_eucledian.order, cluster)) 

## map the clusteres
sync.overlap <- matchClustersToPhase(sync.hc_dtw.df, BD.markers.sig)
sync.overlap$cluster <- as.numeric(sync.overlap$cluster)

sync.phase.match <- sync.overlap %>% group_by(cluster) %>% summarize(phase = markers[which.max(percent)])
sync.tc.mus.scale <- left_join(sync.tc.mus.scale, sync.phase.match, by = 'cluster')

sync.tc.mus.scale$phase <- factor(sync.tc.mus.scale$phase, levels = c('G1', 'S', 'M', 'C'))

p1 <- ggplot(sync.tc.mus.scale, aes(x = t, y = reorder_within(GeneID, -peak.ord, phase), fill = y)) + 
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

plot(p1)  


ggsave(filename="../Output/compScBdTgPb/figs/bd_curve_cluster_heatmap_sync.pdf",
       plot=p1,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

sc.overlap <- matchClustersToPhase(sc.hc_dtw.df, BD.markers.sig)
sc.overlap$cluster <- as.numeric(sc.overlap$cluster)

sc.phase.match <- sc.overlap %>% group_by(cluster) %>% summarize(phase = markers[which.max(percent)])
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.phase.match, by = 'cluster')

sc.tc.mus.scale$phase <- factor(sc.tc.mus.scale$phase, levels = c('G1', 'S', 'M', 'C'))

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


ggsave(filename="../Output/compScBdTgPb/figs/bd_curve_cluster_heatmap_sc.pdf",
       plot=p2,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

p3 <- ggplot(sc.overlap, aes(x = cluster, y = markers, fill = percent)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("markers") + xlab("clusters") + theme_bw() + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(1000)) +
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


saveRDS(sc.tc.mus, '../Input/compScBdTgPb/RData/bd_sc_tc_mus.RData')
saveRDS(sync.tc.mus, '../Input/compScBdTgPb/RData/bd_sync_tc_mus.RData')
saveRDS(sc.tc.mus.scale, '../Input/compScBdTgPb/RData/bd_sc_tc_mus_scale.RData')
saveRDS(sync.tc.mus.scale, '../Input/compScBdTgPb/RData/bd_sync_tc_mus_scale.RData')


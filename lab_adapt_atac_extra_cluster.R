library(dtwclust)
library(bigmemory)
library(doParallel)
library(tidyverse)

library(tidytext)
library(gridExtra)

library(grid)
#library(fda)
#library(sme)
library(openxlsx)
#library(tvReg)

source("./util_funcs.R")

plot.trend <- function(gene.id, gene.name, extra.tc.logCPM){
  
  tmp <- extra.tc.logCPM %>% dplyr::filter(GeneID == gene.id) %>%
    transmute(GeneID = GeneID, x = Time, y = expr, trend.pval = trend.pval, trend.fdr = trend.fdr)
  
  p <- ggplot(tmp, aes(x = x, y = y)) + 
    geom_point(aes(x = x, y = y, color = 'red')) + 
    geom_smooth(aes(x = x, y = y), method='lm', color = 'black') + theme_bw() + 
    ylab('Expr') + xlab('Passage') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste(gene.name, gene.id, sep = ':')) +
    theme(
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) 
  
  
  return(p)
  
}



extra.tc.logCPM <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_atac_logCPM.RData')

gene.id <- 'TGGT1_268850'
gene.name <- 'ENO2'

gene.id <- "TGGT1_227290"
gene.name <- 'MORC'

gene.id <- 'TGGT1_254555'
gene.name <- 'GCN5-A'
p <- plot.trend(gene.id, gene.name, extra.tc.logCPM)

plot(p)

genes <- unique(extra.tc.logCPM$GeneID)


extra.atac.spline.fits <- mclapply(1:length(genes), function(i){
  tmp <- extra.tc.logCPM %>% dplyr::filter(GeneID == genes[i]) %>%
    transmute(GeneID = GeneID, x = Time, y = expr, trend.pval = trend.pval, trend.fdr = trend.fdr)
  extra.atac.sp <- smooth.spline(tmp$x, tmp$y)
  extra.atac.sp <- predict(extra.atac.sp, seq(11, 210,length.out = 10)) 
  mu <- data.frame(x = extra.atac.sp$x, y = extra.atac.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  fit <- lm(mu$y~mu$x)
  res <- summary(fit)
  mu$adj.r.squared <- res$adj.r.squared
  mu$trend.pval <- tmp$trend.pval[1]
  mu$trend.fdr <- tmp$trend.fdr[1]
  
  return(mu)
}, mc.cores = num.cores)


extra.atac.spline.fits <- bind_rows(extra.atac.spline.fits)

extra.atac.spline.fits.filt <- extra.atac.spline.fits %>% dplyr::filter(adj.r.squared > 0.5 & trend.fdr < 0.35)

extra.atac.dtw.wide <- extra.atac.spline.fits.filt %>% 
  pivot_wider(-c(adj.r.squared, trend.fdr, trend.pval), names_from = 'GeneID', values_from = 'y') %>%
  mutate_at(vars(matches('TGGT1')), scale) %>%
  as.data.frame()

## Generate the clusters
num.clust <- 2L

extra.atac.hc_dtws <- dtwClustCurves(extra.atac.dtw.wide[2:ncol(extra.atac.dtw.wide)], nclust = num.clust)

saveRDS(extra.atac.hc_dtws, '../Input/compScBdTgPb/LabAdaptationRNA/extra.atac.hc_dtws_all.rds')

## Run from here
extra.atac.hc_dtws <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra.atac.hc_dtws.rds')

plot(extra.atac.hc_dtws, type = 'sc')
plot(extra.atac.hc_dtws, type = "series", clus = 1L)
plot(extra.atac.hc_dtws, type = "centroids", clus = 1L)


## plot a few curves
## Plot trends

for(i in 1:20){
  tmp <- extra.atac.spline.fits %>% dplyr::filter(GeneID == sig.genes[i])
  plot(tmp$x, tmp$y, type = 'l', lwd = 2, col = 'red')
  Sys.sleep(0.8)
}



atac.clust.info <- data.frame(GeneID = colnames(extra.atac.dtw.wide)[2:ncol(extra.atac.dtw.wide)], 
                             cluster = cutree(extra.atac.hc_dtws, k = 2))

gene.groups <- atac.clust.info %>% group_by(cluster) %>% summarise(genes = list(GeneID))


extra.atac.dtw.long <- extra.atac.dtw.wide %>% pivot_longer(-x, names_to = 'GeneID', values_to = 'y')
extra.atac.dtw.long <- left_join(extra.atac.dtw.long, atac.clust.info, by = 'GeneID')

fits <- lapply(1:nrow(gene.groups), function(i){
  #tmp <- extra.atac.spline.fits[extra.atac.spline.fits$GeneID %in% unlist(gene.groups$genes[i]), ]
  tmp <- extra.atac.dtw.long[extra.atac.dtw.long$GeneID %in% unlist(gene.groups$genes[i]), ]
  fit <- lm(tmp$y~tmp$x)
  summary(fit)
})

atac.r2 <- data.frame(R2 = unlist(lapply(fits, function(x) x$adj.r.squared)), 
                     slope = unlist(lapply(fits, function(x) x$coefficients[2])),
                     cluster = 1:nrow(gene.groups))

extra.atac.dtw.long <- left_join(extra.atac.dtw.long, atac.r2, by = 'cluster')
extra.atac.dtw.long$is.trending <- ifelse(extra.atac.dtw.long$R2 > 0.4, 'yes', 'no')
extra.atac.dtw.long$trending <- ifelse(extra.atac.dtw.long$slope > 0, 'up', 'down')
extra.atac.dtw.long.filt <- extra.atac.dtw.long %>% dplyr::filter(is.trending == 'yes')
extra.atac.dtw.long.filt$trending <- factor(extra.atac.dtw.long.filt$trending, levels = c('up', 'down'))

p <- ggplot(extra.atac.dtw.long.filt, aes(x = x, y = y, group = GeneID)) + 
  geom_line(aes(x = x, y = y, color = trending), alpha = 0.4) + 
  geom_smooth(aes(x = x, y = y, group = NA), method='lm', color = 'black') + 
  theme_bw(base_size = 14) +
  ylab('Access') + xlab('Passage') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  facet_grid(trending~.) + 
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


plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/lab_adapt_trending_atac_clusters.pdf", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


saveRDS(extra.atac.dtw.long, "../Input/compScBdTgPb/RData/extra_atac_dtw_trending_clusters.rds")

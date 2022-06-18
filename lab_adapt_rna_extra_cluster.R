library(bigmemory)
library(dtwclust)
library(doParallel)
library(tidyverse)
library(tidytext)
library(gridExtra)
library(grid)
library(fda)
library(sme)
library(openxlsx)
library(tvReg)

source("./util_funcs.R")

## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


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


extra.tc.logCPM <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_logCPM.RData')

gene.id <- 'TGGT1_268850'
gene.name <- 'ENO2'

gene.id <- 'TGGT1_268860'
gene.name <- 'ENO1'

gene.id <- "TGGT1_227290"
gene.name <- 'MORC'

gene.id <- 'TGGT1_254555'
gene.name <- 'GCN5-A'

p <- plot.trend(gene.id, gene.name, extra.tc.logCPM)

plot(p)

genes <- unique(extra.tc.logCPM$GeneID)

extra.rna.spline.fits <- mclapply(1:length(genes), function(i){
  tmp <- extra.tc.logCPM %>% dplyr::filter(GeneID == genes[i]) %>%
    transmute(GeneID = GeneID, x = Time, y = expr, trend.pval = trend.pval, trend.fdr = trend.fdr)
  extra.rna.sp <- smooth.spline(tmp$x, tmp$y)
  extra.rna.sp <- predict(extra.rna.sp, seq(11, 210,length.out = 10)) 
  mu <- data.frame(x = extra.rna.sp$x, y = extra.rna.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  fit <- lm(mu$y~mu$x)
  res <- summary(fit)
  mu$adj.r.squared <- res$adj.r.squared
  mu$trend.pval <- tmp$trend.pval[1]
  mu$trend.fdr <- tmp$trend.fdr[1]
  
  return(mu)
}, mc.cores = num.cores)

extra.rna.spline.fits <- bind_rows(extra.rna.spline.fits)



extra.rna.spline.fits.filt <- extra.rna.spline.fits %>% dplyr::filter(adj.r.squared > 0.5 & trend.fdr < 0.2)


extra.rna.dtw.wide <- extra.rna.spline.fits.filt %>% 
  pivot_wider(-c(adj.r.squared, trend.fdr, trend.pval), names_from = 'GeneID', values_from = 'y') %>%
  mutate_at(vars(matches('TGGT1')), scale) %>%
  as.data.frame()



## Generate the clusters
num.clust <- 2L


extra.rna.hc_dtws <- dtwClustCurves(extra.rna.dtw.wide[2:ncol(extra.rna.dtw.wide)], nclust = num.clust)

saveRDS(extra.rna.hc_dtws, '../Input/compScBdTgPb/LabAdaptationRNA/extra.rna.hc_dtws.rds')

## Run from here
extra.rna.hc_dtws <- readRDS('../Input/compScBdTgPb/LabAdaptationRNA/extra.rna.hc_dtws.rds')

plot(extra.rna.hc_dtws, type = 'sc')
plot(extra.rna.hc_dtws, type = "series", clus = 1L)
plot(extra.rna.hc_dtws, type = "centroids", clus = 1L)



rna.clust.info <- data.frame(GeneID = colnames(extra.rna.dtw.wide)[2:ncol(extra.rna.dtw.wide)], 
                             cluster = cutree(extra.rna.hc_dtws, k = 4))

gene.groups <- rna.clust.info %>% group_by(cluster) %>% summarise(genes = list(GeneID))


extra.rna.dtw.long <- extra.rna.dtw.wide %>% pivot_longer(-x, names_to = 'GeneID', values_to = 'y')
extra.rna.dtw.long <- left_join(extra.rna.dtw.long, rna.clust.info, by = 'GeneID')

fits <- lapply(1:nrow(gene.groups), function(i){
  #tmp <- extra.rna.spline.fits[extra.rna.spline.fits$GeneID %in% unlist(gene.groups$genes[i]), ]
  tmp <- extra.rna.dtw.long[extra.rna.dtw.long$GeneID %in% unlist(gene.groups$genes[i]), ]
  fit <- lm(tmp$y~tmp$x)
  summary(fit)
})

rna.r2 <- data.frame(R2 = unlist(lapply(fits, function(x) x$adj.r.squared)), 
                     slope = unlist(lapply(fits, function(x) x$coefficients[2])),
                     cluster = 1:nrow(gene.groups))

extra.rna.dtw.long <- left_join(extra.rna.dtw.long, rna.r2, by = 'cluster')
extra.rna.dtw.long$is.trending <- ifelse(extra.rna.dtw.long$R2 > 0.4, 'yes', 'no')
extra.rna.dtw.long$trending <- ifelse(extra.rna.dtw.long$slope > 0, 'up', 'down')
extra.rna.dtw.long.filt <- extra.rna.dtw.long %>% dplyr::filter(is.trending == 'yes')
extra.rna.dtw.long.filt$trending <- factor(extra.rna.dtw.long.filt$trending, levels = c('up', 'down'))


p <- ggplot(extra.rna.dtw.long.filt, aes(x = x, y = y, group = GeneID)) + 
  geom_line(aes(x = x, y = y, color = trending), alpha = 0.4) + 
  geom_smooth(aes(x = x, y = y, group = NA), method='lm', color = 'black') + 
  theme_bw(base_size = 14) +
  ylab('Expr') + xlab('Passage') +
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

ggsave(filename="../Output/compScBdTgPb/figs/lab_adapt_trending_rna_clusters.pdf", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


saveRDS(extra.rna.dtw.long, "../Input/compScBdTgPb/RData/extra_rna_dtw_trending_clusters.rds")

extra.rna.dtw.long <- readRDS("../Input/compScBdTgPb/RData/extra_rna_dtw_trending_clusters.rds")

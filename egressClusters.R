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
library(dtwclust)
library(doParallel)




source('./util_funcs.R')

bd.prod.desc <- read.csv('../Input/compScBdTgPb/genes/BDvi_Prod_desc.csv')
bd.prod.desc <- bd.prod.desc %>% transmute(GeneID = Gene.ID, ProductDescription = Product.Description)

bd.prod.desc[grep('calcium', bd.prod.desc$ProductDescription),]
CDPK4.id <- "Bdiv_024410"
bd.prod.desc[grep('GMP',bd.prod.desc$ProductDescription),]
PKG.id <- "Bdiv_020500"

egress_markers <- read.xlsx('../Input/compScBdTgPb/gene_function/sync_egress_kz.xlsx', sheet = 1)
egress_markers.putative <- read.xlsx('../Input/BdivCellCycle/gene_function/sync_egress_kz.xlsx', sheet = 2)

bd.sc.tc.mus.scale <- readRDS('../Input/compScBdTgPb/RData/bd_sc_tc_mus_scale.RData')
## Generate the clusters
num.clust <- 4L

bd.sc.tc.mus.scale$phase <- factor(sc.tc.mus.scale$phase, levels = c('G1', 'S', 'M', 'C'))
bd.gene.ord <- bd.sc.tc.mus.scale %>% dplyr::select(GeneID, peak.ord) %>% distinct()
bd.sc.tc.mus.scale$GeneID <- factor(bd.sc.tc.mus.scale$GeneID, 
                                    levels = bd.gene.ord$GeneID[order(-bd.gene.ord$peak.ord)])
p <- ggplot(bd.sc.tc.mus.scale, aes(x = t, y = GeneID, fill = y)) + 
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
  scale_y_discrete(breaks = gsub('-', '_', egress_markers$GeneID), labels = egress_markers$Bd_name, 
                   guide=guide_axis(n.dodge=3)) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=5, face="bold", color = 'black')
  ) 

plot(p)  


ggsave(filename="../Output/compScBdTgPb/figs/bd_curve_cluster_heatmap_sc_egress.pdf",
       plot=p,
       width = 6, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)




## Recluster the egress markers
bd.sc.tc.mus.scale$GeneID <- as.character(bd.sc.tc.mus.scale$GeneID)

bd.sc.egress.markers.tc.mus <- bd.sc.tc.mus.scale %>% dplyr::filter(GeneID %in% gsub('-', '_', egress_markers$GeneID))

bd.sc.egress.markers.tc.mus <- bd.sc.egress.markers.tc.mus %>% dplyr::select(GeneID, t, y) %>% 
  pivot_wider(names_from = GeneID, values_from = y)  %>% as.data.frame()

bd.sc.egress.markers.hc_dtw <- dtwClustCurves(bd.sc.egress.markers.tc.mus[,2:ncol(bd.sc.egress.markers.tc.mus)], nclust = 4L)

plot(bd.sc.egress.markers.hc_dtw, type = 'sc')

bd.sc.egress.markers.tc.mus.long <- bd.sc.egress.markers.tc.mus %>% 
  pivot_longer(cols = -t, names_to = 'GeneID', values_to = 'y')

clust.info <- data.frame(GeneID = colnames(bd.sc.egress.markers.tc.mus[,2:ncol(bd.sc.egress.markers.tc.mus)]), 
                         cluster = cutree(bd.sc.egress.markers.hc_dtw, k = 4))
bd.sc.egress.markers.tc.mus.long <- inner_join(bd.sc.egress.markers.tc.mus.long, clust.info, by = 'GeneID')
bd.sc.egress.markers.tc.mus.long$Name <- egress_markers$Bd_name[match(bd.sc.egress.markers.tc.mus.long$GeneID, 
                                                                   gsub('-', '_', egress_markers$GeneID))]


p <- lapply(unique(bd.sc.egress.markers.tc.mus.long$cluster), function(i){
  tmp <- bd.sc.egress.markers.tc.mus.long %>% dplyr::filter(cluster == i)
  ggplot(tmp, aes(x = t, y = y, color = Name)) + 
    geom_line() + theme_bw()
})

p2 <- do.call(grid.arrange, c(p, nrow=2))

ggsave(filename="../Output/compScBdTgPb/figs/bd_sc_egress_cluster_curves.pdf",
       plot=p2,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)



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



saveRDS(bd.sc.egress.markers.tc.mus.long, '../Input/compScBdTgPb/RData/bd_sc_egress._markers_tc_mus_long.RData')

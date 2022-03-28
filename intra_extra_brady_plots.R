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
library(sctransform)
library(cowplot)
library(patchwork)
library(doParallel)
library(ggVennDiagram)
library(tidytext)


source('./util_funcs.R')

## Extracting the data for control over PCA directions
getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp2, phase = S.O@meta.data$phase)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

S.Os <- readRDS('../Input/compScBdTgPb/RData/S.O.intra_extra_pruWT_pruPH_lables_list.RData')
S.O.integrated <- readRDS('../Input/compScBdTgPb/RData/S.O.intra.extra.pruWT.pruPH.integrated.rds')
S.O.no.anchored <- readRDS('../Input/compScBdTgPb/RData/S.O.intra.extra.pruWT.pruPH.no_anchor.rds')


TF.info <- read.xlsx('../Input/compScBdTgPb/genes/TF_Info_Updated_kz.xlsx')
genes <- gsub("_", "-", TF.info$GeneID)
labs <- gsub("_", "-", TF.info$Ap2Name)


### Adding better spp names
S.O.integrated@meta.data$spp2 <- ifelse(S.O.integrated@meta.data$spp == 'intra', 'RH.intra',
                                        ifelse(S.O.integrated@meta.data$spp == 'extra', 'RH.extra',
                                               ifelse(S.O.integrated@meta.data$spp == 'pru.WT', 'Pru.intra',
                                                      'Pru.brady')))
S.O.integrated@meta.data$spp2 <- factor(S.O.integrated@meta.data$spp2, 
                                        levels = c('RH.intra', 'RH.extra', 'Pru.intra', 'Pru.brady'))
S.O.integrated$phase.spp2 <- paste(S.O.integrated@meta.data$spp2, S.O.integrated@meta.data$phase, sep = "_")
Idents(S.O.integrated) <- "phase.spp2"


S.O.no.anchored@meta.data$spp2 <- ifelse(S.O.no.anchored@meta.data$spp == 'intra', 'RH.intra',
                                         ifelse(S.O.no.anchored@meta.data$spp == 'extra', 'RH.extra',
                                                ifelse(S.O.no.anchored@meta.data$spp == 'pru.WT', 'Pru.intra',
                                                       'Pru.brady')))
S.O.no.anchored@meta.data$spp2 <- factor(S.O.no.anchored@meta.data$spp2, 
                                         levels = c('RH.intra', 'RH.extra', 'Pru.intra', 'Pru.brady'))
Idents(S.O.no.anchored) <- 'spp2'


## Violin plots of all known regulator
p <- VlnPlot(S.O.no.anchored, features = genes)

out.dir <- "../Output/compScBdTgPb/figs/expression_plots/"
p <- lapply(1:length(genes), function(i){
  pp <- p[[i]] + 
    ylab('Expr') + xlab('') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    #facet_wrap(spp~.) + 
    ggtitle(paste(labs[i], gsub('-', '_', genes[i]),sep = ':')) +
    theme(
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) + 
    theme(#legend.position = c(0.1, 0.25),
      legend.position = 'None',
      legend.title = element_text(colour="black", size=6, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=6, 
                                 face="bold"))
  
  f.name <- paste(paste(out.dir, paste(labs[i], gsub('-', '_', genes[i]),sep = '_'), sep = ''), '.pdf', sep = '')
  ggsave(filename=f.name,
         plot=pp,
         width = 5, height = 4,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
  )
  
  pp
})


## UMAP plots
pcaMataData.alldata <- getPcaMetaData(S.O.no.anchored)

pcaMataData.alldata$phase <- factor(pcaMataData.alldata$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))

p2  <- ggplot(pcaMataData.alldata, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = spp,
    color = spp, 
  ), #color = 'blue', 
  alpha = 0.9,
  shape=21, size = 0.3)+ 
  #scale_color_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  #scale_fill_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('UMAP2') + xlab('UMAP1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = c(0.13, 0.85),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))


plot(p2)

ggsave(filename='../Output/compScBdTgPb/figs/projection_plots/umap_not_anchored.pdf',
       plot=p2,
       width = 6.5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


  pcaMataData.S.O.integrated<- getPcaMetaData(S.O.integrated)
pcaMataData.S.O.integrated$phase <- factor(pcaMataData.S.O.integrated$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))


p5  <- ggplot(pcaMataData.S.O.integrated, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase, 
  ), #color = 'blue', 
  alpha = 0.9,
  shape=21, size = 0.3)+ 
  #scale_color_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  #scale_fill_manual(values = c("intra" = "firebrick","extra" ="darkorchid3", 'pru.pH' = 'darkslateblue')) +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('UMAP2') + xlab('UMAP1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  facet_wrap(spp~.) + 
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


plot(p5)


ggsave(filename='../Output/compScBdTgPb/figs/projection_plots/umap_anchored_split_phase.pdf',
       plot=p5,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## markers
Global.Markers.sig <- read.xlsx('../Output/compScBdTgPb/tables/Tachy_Brady_Extra_global_markers.xlsx')

Global.Markers.stat <- Global.Markers.sig %>% group_by(stage) %>% summarise(num.deg = n())
Global.Markers.stat$stage <- factor(Global.Markers.stat$stage,
                                    levels = c('Tachy', 'Extra', 'Brady'))
# bar plot 
p <-  ggplot(Global.Markers.stat, aes(x=stage, y=num.deg)) +
  geom_bar(stat="identity", fill = "steelblue")+
  geom_text(aes(label=num.deg), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Intra (RH/Pru)", "Extra (RH)", "Brady (Pru)")) +
  theme(axis.text.x = element_text(face="bold", size=14, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=14, angle=0, color = 'black'),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/Tachy_Extra_Brady_unique_markers.pdf", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


Brady.vs.Tachy.markers.sig <- read.xlsx('../Output/compScBdTgPb/tables/Brady_vs_Tachy_markers.xlsx')

Brady.vs.Tachy.markers.sig$stage <- ifelse(Brady.vs.Tachy.markers.sig$avg_log2FC > 0, 'Brady', 'Tachy')
Brady.vs.Tachy.markers.stat <- Brady.vs.Tachy.markers.sig %>% group_by(stage) %>% summarise(num.deg = n())
Brady.vs.Tachy.markers.stat$stage <- factor(Brady.vs.Tachy.markers.stat$stage,
                                    levels = c('Tachy', 'Brady'))
# bar plot 
p <-  ggplot(Brady.vs.Tachy.markers.stat, aes(x=stage, y=num.deg)) +
  geom_bar(stat="identity", fill = "steelblue")+
  geom_text(aes(label=num.deg), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Intra (RH/Pru)",  "Brady (Pru)")) +
  theme(axis.text.x = element_text(face="bold", size=14, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=14, angle=0, color = 'black'),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/Tachy_Brady_total_DEGs.pdf", 
       plot=p,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

Extra.vs.Intra.RH.markers.sig <- read.xlsx('../Output/compScBdTgPb/tables/Extra_vs_Intra_RH_markers.xlsx')
Extra.vs.Intra.RH.markers.sig$stage <- ifelse(Extra.vs.Intra.RH.markers.sig$avg_log2FC > 0, 'Extra', 'Intra')
Extra.vs.Intra.RH.markers.stat <- Extra.vs.Intra.RH.markers.sig %>% group_by(stage) %>% summarise(num.deg = n())
Extra.vs.Intra.RH.markers.stat$stage <- factor(Extra.vs.Intra.RH.markers.stat$stage,
                                            levels = c('Intra', 'Extra'))
# bar plot 
p <-  ggplot(Extra.vs.Intra.RH.markers.stat, aes(x=stage, y=num.deg)) +
  geom_bar(stat="identity", fill = "steelblue")+
  geom_text(aes(label=num.deg), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Intra (RH)",  "Extra (Pru)")) +
  theme(axis.text.x = element_text(face="bold", size=14, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=14, angle=0, color = 'black'),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/Intra_Extra_total_DEGs.pdf", 
       plot=p,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Shared markers
shared.Brady.Extra <- read.xlsx('../Output/compScBdTgPb/tables/shared_brady_extra_markers.xlsx')
## Using Intra as base
brady.up <- Brady.vs.Tachy.markers.sig %>% dplyr::filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
brady.down <- Brady.vs.Tachy.markers.sig %>% dplyr::filter(avg_log2FC < 0) %>% arrange(avg_log2FC)
Extra.up <- Extra.vs.Intra.RH.markers.sig %>% dplyr::filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
Extra.down <- Extra.vs.Intra.RH.markers.sig %>% dplyr::filter(avg_log2FC < 0) %>% arrange(avg_log2FC)

Brady.Extra <- list(Brady.up = brady.up$GeneID, Brady.down = brady.down$GeneID, 
                    Extra.up = Extra.up$GeneID, Extra.down = Extra.down$GeneID)

g <- ggVennDiagram(Brady.Extra, label_size = 5, set_size = 6,label = "count", label_alpha = 0,label_color = 'white',
                   category.names = c("Brady (Up)", "Brady (Down)", "Extra (Up)", "Extra (Down)"),
                   set_color = c("brady.up" = "firebrick","brady.down" ="darkorchid3", 'Extra.up' = 'darkslateblue', 
                                 'Extra.down' = 'darkolivegreen4')) +  
  scale_color_manual(values = c("firebrick","darkorchid3", 'darkslateblue', 'darkolivegreen4')) + 
  scale_x_continuous(expand = expansion(mult = .15)) + theme(legend.position = "none")

plot(g)


ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/Brady_Extra_shared_degs.pdf", 
       plot=g,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Expression heatmaps
exprs <- as.matrix(S.O.integrated.sig[["RNA"]]@data) %>% data.frame()
exprs$gene <- rownames(exprs)
exprs.sig <- exprs %>% dplyr::filter(gene %in% gsub('_', '-', TF.info$GeneID))

exprs.sig.ave <- exprs.sig %>% transmute(GeneID = gsub('-', '_', gene),
                                         intra = rowMeans(across(matches("intra"))),
                                         extra = rowMeans(across(contains("extra"))),
                                         brady = rowMeans(across(contains("pru.pH"))))

exprs.sig.ave <- exprs.sig.ave %>% mutate(stage = names(across(matches("intra|extra|brady")))[max.col(across(matches("intra|extra|brady")))])

heatmap(as.matrix(exprs.sig.ave[,-c(1)]), scale = "col", Rowv = NA, Colv = NA)

# hc_eucledian <- hclust(dist((as.matrix(exprs.sig.ave[,-c(1,2)] ))), method = "ward.D")
# hc_eucledian.df <- data.frame(GeneID = gsub('-', '_', hc_eucledian$labels),
#                               gene.ord = hc_eucledian$order, gene.cluster = cutree(hc_eucledian,k = 3))

eexprs.sig.ave.long <- exprs.sig.ave %>% pivot_longer(-c(GeneID,gene.ord),
                                                     names_to = 'stage', values_to = 'ave_expr')
exprs.sig.ave.long$GeneID <- factor(exprs.sig.ave.long$GeneID, levels = exprs.sig.ave$GeneID[exprs.sig.ave$gene.ord])
p <- ggplot(exprs.sig.ave.long, aes(x = stage, y = reorder(GeneID, gene.ord),  fill = ave_expr)) +
  geom_tile() +
  scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("stage") +
  # facet_grid(stage.marker~stage, scales = 'free', space = 'free')+
  # theme(strip.text = element_text(size = 18, face = 'bold')) + 
  # scale_fill_gradientn(colours = viridis::inferno(10)) +
  # theme(
  #   axis.text.x = element_blank(),
  #   #axis.text.y = element_text(size = 18),
  #   axis.text.y = element_blank(),
  #   axis.ticks = element_blank(),
  #   axis.title.x = element_text(size = 18, face = "bold"),
  #   axis.title.y = element_text(size = 18, face = "bold"),
  #   legend.position = "none")+
  theme(strip.background=element_rect(fill='white', color = 'black'))+
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(strip.text.y=element_text(angle=0, hjust=0.5,vjust=0.5, face = 'bold'))+
  #ggtitle(titles)+
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  )
p


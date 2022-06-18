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
library(dittoSeq)



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

S.O.intra.extra <- readRDS('../Input/compScBdTgPb/RData/S.O.intra.extra.labs_list.rds')

S.O.intra <- S.O.intra.extra$intra
S.O.intra@meta.data$spp2 <- S.O.intra@meta.data$spp
intra.pca <- getPcaMetaData(S.O.intra)

S.O.extra <- S.O.intra.extra$extra
S.O.extra@meta.data$spp2 <- S.O.extra@meta.data$spp
extra.pca <- getPcaMetaData(S.O.extra)

p1  <- ggplot(extra.pca, aes(x= -UMAP_1,y=-UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  
  theme_bw(base_size = 14) +
  #theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  #ylab('UMAP2') + xlab('UMAP1') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  # facet_grid(.~spp) + 
  # theme(
  #   strip.text.x = element_text(
  #     size = 18, color = "black", face = "bold.italic"
  #   )
  # ) + 
  theme(legend.position="none") + 
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 0),
    axis.title.y = element_text(size=18, face="bold")
  )


plot(p1)




ggsave(filename='../Output/compScBdTgPb/figs/umap_extra.pdf',
       plot=p1,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)





S.Os <- readRDS('../Input/compScBdTgPb/RData/S.O.intra_extra_pruWT_pruPH_lables_list.RData')
S.O.integrated <- readRDS('../Input/compScBdTgPb/RData/S.O.intra.extra.pruWT.pruPH.integrated.rds')
S.O.no.anchored <- readRDS('../Input/compScBdTgPb/RData/S.O.intra.extra.pruWT.pruPH.no_anchor.rds')



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



integratred.pca <- getPcaMetaData(S.O.integrated)
p1  <- ggplot(integratred.pca, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  
  theme_bw(base_size = 14) +
  #theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  #ylab('UMAP2') + xlab('UMAP1') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  # facet_grid(.~spp) + 
  # theme(
  #   strip.text.x = element_text(
  #     size = 18, color = "black", face = "bold.italic"
  #   )
  # ) + 
  theme(legend.position = "None",
        #legend.position = c(0.84, 0.14),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))+ 
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 0),
    axis.title.y = element_text(size=18, face="bold")
  )


plot(p1)




ggsave(filename='../Output/compScBdTgPb/figs/umap_integrated_lourido.pdf',
       plot=p1,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



S.Os.boothroyd <- readRDS('../Input/compScBdTgPb/RData/S.O.intra_extra_pru_d0_d7_lables_list_boothroyd.RData')
S.O.integrated.boothroyd <- readRDS('../Input/compScBdTgPb/RData/S.O.intra.extra.pru_d0_d7.integrated_boothroyd.rds')

integratred.pca.boothroyd <- getPcaMetaData(S.O.integrated.boothroyd)
p1  <- ggplot(integratred.pca.boothroyd, aes(x= -UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  
  theme_bw(base_size = 14) +
  #theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  #ylab('UMAP2') + xlab('UMAP1') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  # facet_grid(.~spp) + 
  # theme(
  #   strip.text.x = element_text(
  #     size = 18, color = "black", face = "bold.italic"
  #   )
  # ) + 
  theme(legend.position = "None",
        #legend.position = c(0.84, 0.14),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=3)))+ 
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 0),
    axis.title.y = element_text(size=18, face="bold")
  )


plot(p1)




ggsave(filename='../Output/compScBdTgPb/figs/umap_integrated_boothroyd.pdf',
       plot=p1,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## Violin plots of all known regulator
TF.info <- read.xlsx('../Input/compScBdTgPb/genes/TF_Info_Updated_kz.xlsx')
genes <- gsub("_", "-", TF.info$GeneID)
labs <- gsub("_", "-", TF.info$Ap2Name)


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



## Violin plot for stage conversion players (from Sokol, Boyle review)
Stage.info <- read.xlsx('../Input/compScBdTgPb/gene_function/Brady_Stage_Conversion_Sokol_Boyle_review.xlsx')
genes <- gsub("_", "-", Stage.info$GeneID)
labs <- gsub("_", "-", Stage.info$GeneName)

DefaultAssay(S.O.integrated) <- 'RNA'
p <- VlnPlot(S.O.no.anchored, features = genes)

out.dir <- "../Output/compScBdTgPb/figs/stage_conversion_plots/"
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



## Dot plots
DefaultAssay(S.O.integrated) <- 'RNA'

final.ap2s <- c('AP2Ib-1', 
                'AP2IX-4', 
                'AP2IX-9', 
                'AP2VIII-7', 
                'AP2IX-3', 
                'AP2VIIa-6', 
                'TFIIS', 
                'AP2IX-10')

final.ap2.ids <- gsub('_', '-', TF.info$GeneID[match(final.ap2s, TF.info$Ap2Name)])

DotPlot(
  S.O.integrated, features = final.ap2.ids,
  group.by = "spp2") + coord_flip() + theme_light() + 
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) + 
  scale_x_discrete(labels = final.ap2s) + 
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1, color = 'black'),
    axis.title.y = element_text(size=14, face="bold", color = 'black')
  ) + 
  theme(#legend.position = c(0.1, 0.25),
    legend.position = 'None',
    legend.title = element_text(colour="black", size=6, 
                                face="bold"),
    legend.text = element_text(colour="black", size=6, 
                               face="bold"))
out.dir <- '../Output/compScBdTgPb/figs/final_ap2s/'
Idents(S.O.integrated) <- 'spp2'
p <- VlnPlot(S.O.integrated, features = final.ap2.ids)
p <- lapply(1:length(final.ap2.ids), function(i){
  pp <- p[[i]] + 
    ylab('') + xlab('') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    #facet_wrap(spp~.) + 
    ggtitle(paste(final.ap2s[i], gsub('-', '_', final.ap2.ids[i]),sep = ':')) +
    theme(
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(size=12, face="bold", hjust = 0.5),
      axis.title.y = element_text(size=12, face="bold")
    ) + 
    theme(#legend.position = c(0.1, 0.25),
      legend.position = 'None',
      legend.title = element_text(colour="black", size=6, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=6, 
                                 face="bold"))
  
  f.name <- paste(paste(out.dir, paste(final.ap2s[i], gsub('-', '_', final.ap2.ids[i]),sep = '_'), sep = ''), '.pdf', sep = '')
  ggsave(filename=f.name,
         plot=pp,
         width = 4, height = 3,
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
  scale_x_discrete(labels=c("Intra (RH)",  "Extra (RH)")) +
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



## Up or Down Barplot
X1 <- Brady.vs.Tachy.markers.sig %>% summarise(stage = 'Brady', total = n())
X2 <- Extra.vs.Intra.RH.markers.sig %>% summarise(stage = 'Extra', total = n())
X3 <- shared.Brady.Extra %>% summarise(stage = 'Shared', total = n())
X <- rbind(X1, X2, X3)

p <-  ggplot(X, aes(x=stage, y=total, fill = stage)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=total), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('DEG') + xlab('') +
  scale_x_discrete(labels=c("Pru (pH 8.4)",  "e.c", 'shared')) +
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
        axis.ticks = element_blank()) + 
  theme(legend.position = "None",
        #legend.position = c(0.88, 0.17),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=2)))



plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/marker_stats/total_up_down_brady_extra_shared.pdf", 
       plot=p,
       width = 4, height = 4, 
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



##### Feature plots
## Other markers
my.genes <- c("TGGT1-208020",
              "TGGT1-240900",
              "TGGT1-280460",
              "TGGT1-306620",
              "TGGT1-214960",
              "TGGT1-203050",
              "TGGT1-309410",
              "TGGT1-264485",
              "TGGT1-269010",
              "TGGT1-215895",
              'TGGT1-200385',
              "TGGT1-259020",
              'TGGT1-233460',
              'TGGT1-311100',
              'TGGT1-209985')

gene.names <- c("AP2Ib-1",
                "AP2VI-2",
                "AP2VIIa-2",
                'AP2IX-9',
                "AP2X-8",
                "AP2VIIa-6",
                "AP2XI-1",
                "AP2IX-3",
                "AP2VIII-7", 
                "AP2IX-10",
                "BFD1",
                "BAG1",
                "SAG1",
                "Zn fingure",
                "cAMP")

i <- 15
p <- FeaturePlot(object = S.O.no.anchored, features = my.genes[i], label = F,repel = T,label.size = 6, 
                 #shape.by  = 'spp',
                 #split.by = 'spp',
                 cols = c("grey", "blue"), reduction = "umap") + 
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
  ggtitle(paste(gene.names[i], gsub('-', '_', my.genes[i]),sep = ':')) +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = c(0.15, 0.8),
        #legend.position = 'None',
        legend.title = element_text(colour="black", size=6, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=6, 
                                   face="bold"))




plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/Zn_Fingure_expr.pdf", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

TGGT1_207600
p <- VlnPlot(S.O.no.anchored, features = 'TGGT1-209985')
plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/Zn_Fingure_expr_violin.pdf", 
       plot=p,
       width = 6, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


##
shared.Brady.Extra <- read.xlsx('../Output/compScBdTgPb/tables/shared_brady_extra_markers_all_genes.xlsx')
shared.Brady.Extra <- left_join(shared.Brady.Extra, prod.desc, by = 'GeneID')

shared.Brady.Extra <- shared.Brady.Extra %>% 
  mutate(up.brady = ifelse(avg_log2FC.Brady > 1, T, F),
         up.extra = ifelse(avg_log2FC.Extra > 1, T, F))


shared.Brady.Extra <- shared.Brady.Extra %>% 
  mutate(color.class = ifelse(avg_log2FC.Brady > log2(1.7)  & avg_log2FC.Extra > log2(1.7) , 'shared', 
                              ifelse(avg_log2FC.Brady > 4 & avg_log2FC.Extra < 0.5, 'brady', 
                                     ifelse(avg_log2FC.Extra >= 2 & avg_log2FC.Brady < 0.6, 'extra','other'))))

shared.Brady.Extra$geneLable <- ''
shared.Brady.Extra <- shared.Brady.Extra %>% 
  mutate(GeneLable = ifelse(!(color.class %in%  c('other')), 'yes', 'no'))

#shared.Brady.Extra$Name <- shared.Brady.Extra$GeneID
shared.Brady.Extra$Name <- paste(shared.Brady.Extra$GeneID, shared.Brady.Extra$ProductDescription, sep = '\n')

shared.Brady.Extra$Name <- gsub('kinase', 'kin.', gsub('serine/threonine prot. phosphatase', 'ser./thre. phosph. ', gsub('Rhoptry kinase fam. prot.', 'Rhop. kinase', gsub('tryptophanyl-tRNA synthetase \\(TrpRS2\\)', 'TrpRS2', 
                                gsub('putative cell-cycle-associated protein kinase PRP4', 'PRP4', 
                                     gsub('hop. kinase fam. prot.', 'Rhop. kinase ', 
                                          gsub('AP2 domain-containing protein', 'AP2 dom.', 
                                               gsub('bradyzoite rhoptry protein ', '', 
                                                    gsub('threonine specific protein phosphatase', 'threonine prot. phosphatase', 
                                                         gsub('cAMP-dependent protein kinase', 'cAMP-dep kinase', 
                                gsub('universal', 'univ.', gsub('phosphatidylinositol 3- and 4-kinase', 'kinase', 
                                gsub('putative cell-cycle-associated protein kinase GSK', 'kinase GSK', 
                                     gsub('bradyzoite antigen ', '', 
                                gsub('putative rhoptry kinase fam. prot. ROP34', 'ROP34', 
                                     gsub('protease inhibitor PI2', 'protease inh. PI2', 
                                          gsub('subtilisin SUB1', 'SUB1', 
                                               gsub('lactate dehydrogenase ', '', 
                                                    gsub('tetratricopeptide repeat-containing protein', 'tetratricopeptide', 
                                                         gsub('putative microneme protein', 'microneme', 
                                                              gsub('AP2 domain transcription factor ', '', 
                                                                   gsub('putative myosin heavy chain', 'GRA46',
                                gsub('zinc finger \\(CCCH type\\) motif-containing protein', 'Zn fing.', 
                                gsub('transporter', 'trans.', gsub('SAG-related sequence ', '', 
                                gsub('family protein', 'fam. prot.', 
                                     gsub('OTU family cysteine protease', 'cysteine protease', 
                                          gsub('hypothetical protein', 'hypo.',
                                gsub('tubulin/FtsZ family, GTPase domain-containing protein', 'tub epsilon',
                                shared.Brady.Extra$Name)))))))))))))))))))))))))))))

shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_311100')] <- 'TGGT1_311100 \n BFD2'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_215210')] <- 'TGGT1_215210 \n F-Box'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_239748')] <- 'TGGT1_239748 \n GPI-GlcNAc transf. \n (fitness!)'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_215970')] <- 'TGGT1_215970 \n hypo. \n excreted'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_247530')] <- 'TGGT1_247530 \n hypo. 1xTM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_203300')] <- 'TGGT1_203300 \n hypo. 2xTM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_217530')] <- 'TGGT1_217530 \n GRA63'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_265320')] <- 'TGGT1_265320 \n splicing factor CWC2'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_311230')] <- 'TGGT1_311230 \n mucin? TM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_273320')] <- 'TGGT1_273320 \n hypo. SP & TM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_321530')] <- 'TGGT1_321530 \n Catheptsin L'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_276170')] <- 'TGGT1_276170 \n PI3,4K'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_312330')] <- 'TGGT1_312330 \n RNA-PolII-B1'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_207160')] <- 'TGGT1_207160 \n SRS49D'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_268860')] <- 'TGGT1_268860 \n ENO1'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_248990')] <- 'TGGT1_248990 \n put. prenylcysteine oxidase'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_272370')] <- 'TGGT1_272370 \n hypo. Sx TM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_223550')] <- 'TGGT1_223550 \n hypo.I 4xTM'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_240090')] <- 'TGGT1_240090 \n ROP34'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_208730')] <- 'TGGT1_208730 \n MIC (putative)'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_225540')] <- 'TGGT1_225540 \n hypo. (TolA?)'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_216140')] <- 'TGGT1_216140 \n TPR domain \n (pub)'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_269690')] <- 'TGGT1_269690 \n hypo. SP'
shared.Brady.Extra$Name[which(shared.Brady.Extra$GeneID == 'TGGT1_253330')] <- 'TGGT1_253330 \n BPK1 (ROP kin)'












shared.Brady.Extra$Name[shared.Brady.Extra$GeneLable == 'yes']
p <- ggplot(shared.Brady.Extra) + geom_point(aes(avg_log2FC.Brady, avg_log2FC.Extra,color=color.class))+ 
  #xlim(-4, 4) + ylim(-3, 3) + 
  labs(color = "") + 
  geom_text_repel(aes(avg_log2FC.Brady, avg_log2FC.Extra), size = 1.7, fontface = "bold",
                  label = ifelse(shared.Brady.Extra$GeneLable == 'yes', as.character(shared.Brady.Extra$Name),""), 
                  box.padding = unit(0.6, "lines"),
                  max.overlaps = 300,
                  segment.angle = 180,
                  nudge_x = 0.5, 
                  nudge_y = 0.5,
                  hjust=0,
                  #nudge_x=0.25, 
                  segment.size = 0.25) + 
  geom_hline(yintercept = log2(1.7), color = 'black', linetype=2) + 
  geom_vline(xintercept = log2(1.7), color = 'black', linetype=2) + 
  xlab("log2 FC Brady over Intra") + ylab("log2 FC Extra over Intra") +
  theme(legend.title=element_blank(),text = element_text(size=20))+ 
  scale_color_manual(values = c("shared" = "red", "brady" = "blue", 'extra' = 'green', 'other' = 'gray'))+
  theme_bw()+
  theme(
    #plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    legend.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.45, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    strip.text = element_text(color = "black", size = 18, face = "bold"),
    strip.background = element_rect(fill = "white"))+
  theme(legend.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.position = c(0.08, 0.25))

p

ggsave(filename="../Output/compScBdTgPb/figs/brady_over_intra_vs_extra_over_intra_fc1_7.pdf", 
       plot=p,
       width = 9, height = 9, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



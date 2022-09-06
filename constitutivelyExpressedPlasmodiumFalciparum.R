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
#library(sctransform)


source('./util_funcs.R')


### P Berghei
input.dir.10x <- "../Input/MalariaCellAtlas/Expression_Matrices/10X/pf10xIDC/"
plasmodium.10x.count.file <- 'pf10xIDC_counts.csv'
plasmodium.10x.pheno.file <- 'pf10xIDC_pheno.csv'


## Reading Plasmodium data 10x pberghei
plasmodium.10x.count <- read.csv(paste(input.dir.10x, plasmodium.10x.count.file, sep = ''))
plasmodium.10x.pheno <- read.csv(paste(input.dir.10x, plasmodium.10x.pheno.file, sep = ''))

genes.10x <- plasmodium.10x.count$X
pf.10x.expr <- plasmodium.10x.count[,-1]
rownames(pf.10x.expr) <- genes.10x

# Figure 3B
plasmodium.10x.pheno <- plasmodium.10x.pheno %>% 
  mutate(cells = case_when(absclust == 0 ~ "TrpE",
                           absclust == 1 ~ "TrpM",
                           absclust == 2 ~ "RngL",
                           absclust == 3 ~ "TrpL",
                           absclust == 4 ~ "SchE",
                           absclust == 5 ~ "SchL",
                           absclust == 6 ~ "RngE",
                           absclust == 7 ~ "SchM"))



plasmodium.10x.pheno$spp <- 'Pf'
plasmodium.10x.pheno$NAME <- paste(plasmodium.10x.pheno$cells, plasmodium.10x.pheno$X, sep = '_')

plasmodium.pheno <- plasmodium.10x.pheno %>% transmute(Sample = X, spp = spp, cells = cells, Name = NAME)

# Set initial Seurat clusters
S.O.pf <- CreateSeuratObject(counts = pf.10x.expr, min.cells = 10, min.features = 100)


VlnPlot(S.O.pf, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.pf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

S.O.pf <- subset(S.O.pf, subset = nFeature_RNA > 150 & nFeature_RNA < 1700 )
S.O.pf <- prep_S.O(S.O.pf)
plasmodium.pheno <- plasmodium.pheno %>% dplyr::filter(Sample %in% colnames(S.O.pf))
rownames(plasmodium.pheno) <- plasmodium.pheno$Sample
S.O.pf <- AddMetaData(S.O.pf, plasmodium.pheno)
Idents(S.O.pf) <- 'cells'

count.data <- S.O.pf@assays$RNA@data
percent.cell.expr <- rowSums(count.data > 0) / ncol(count.data)
ave.expr <- rowSums(count.data)/ncol(count.data)
dd <- data.frame(GeneID = rownames(count.data), percent.cell.expr = percent.cell.expr, ave.expr = ave.expr)
plot(dd$percent.cell.expr, dd$ave.expr, xlim = c(0, 0.8), ylim = c(0, 10))
plot(dd$percent.cell.expr, dd$ave.expr)
q1 = quantile(dd$percent.cell.expr, probs = 0.98)
q2 = quantile(dd$ave.expr, probs = 0.98)

## Top genes
xx <- dd[dd$percent.cell.expr >= q1 & dd$ave.expr >= q2,] %>% arrange(desc(ave.expr))
xx$GeneID <- gsub('-', '_', xx$GeneID)
prod.desc <- read.csv('../Input/compScBdTgPb/genes/Pf_Prod_Desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)
pf.gap <- read.xlsx('../Input/compScBdTgPb/genes/Pf_upstream_intergenic.xlsx')
xx <- left_join(xx, prod.desc, by = 'GeneID')
xx$upstream.intergenic <- pf.gap$upstream[match(xx$GeneID, pf.gap$GeneID)]
xx$strand <- pf.gap$strnd[match(xx$GeneID, pf.gap$GeneID)]
xx <- xx %>% arrange(desc(percent.cell.expr))


write.xlsx(xx, '../Output/compScBdTgPb/tables/pf_high_percent_high_ave_expr.xlsx')

dd$color.class <- 0
dd$color.class[dd$percent.cell.expr >= q1 & dd$ave.expr >= q2] <- 1
dd$color.class <- factor(dd$color.class, levels = c(0,1))

library(ggrepel)
p <- ggplot(dd) + geom_point(aes(percent.cell.expr, ave.expr,color=color.class))+ 
  #xlim(-4, 4) + ylim(-3, 3) + 
  labs(color = "") + 
  geom_text_repel(aes(percent.cell.expr, ave.expr), size = 2, fontface = "bold",
                  label = ifelse(dd$color.class == 1, as.character(dd$GeneID),""), 
                  box.padding = unit(0.6, "lines"),
                  max.overlaps = 300) + 
  # direction = 'both',
  # #segment.angle = 180,
  #nudge_x = 0.5
  # nudge_y = 0.5,
  # hjust=0.25,
  # force = 1,
  # force_pull = 1,
  # segment.size = 0.25) + 
  # 
  # geom_text_repel(aes(percent.cell.expr, ave.expr), size = 2, fontface = "bold",
  #                 label = ifelse(dd$color.class == 1, as.character(dd$GeneID),""), 
#                 box.padding = unit(0.6, "lines"),
#                 max.overlaps = 30,
#                 direction = 'both',
#                 #segment.angle = 180,
#                 nudge_x = 0.5, 
#                 nudge_y = 0.5,
#                 hjust=0.25,
#                 force = 1,
#                 force_pull = 1,
#                 segment.size = 0.25) + 
geom_hline(yintercept = q2, color = 'black', linetype=2) + 
  geom_vline(xintercept = q1, color = 'black', linetype=2) + 
  xlab("percent cell exprssing") + ylab("ave log2(expr)") +
  theme(legend.title=element_blank(),text = element_text(size=20))+ 
  scale_color_manual(values = c("1" = "red", "0" = 'gray'))+
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
        #legend.position = c(0.08, 0.25),
        legend.position = 'none')

p

ggsave(filename="../Output/compScBdTgPb/figs/pf_percent_cell_vs_ave_expr.pdf", 
       plot=p,
       width = 12, height = 12, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



Idents(S.O.pf) <- 'cells'
p <- FeaturePlot(S.O.pf, features = gsub('_', '-', xx$GeneID[1:12]), reduction = 'umap', label = T, label.size = 5)
plot(p)
ggsave(filename="../Output/compScBdTgPb/figs/pf_expr_high_percent_high_ave_expr.pdf", 
       plot=p,
       width = 12, height = 10, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## calculate the overlap
xx.pf <- read.xlsx('../Output/compScBdTgPb/tables/pf_high_percent_high_ave_expr.xlsx')
xx.pb <- read.xlsx('../Output/compScBdTgPb/tables/pb_high_percent_high_ave_expr.xlsx')
pf.pb.orth <- read.xlsx('../Input/compScBdTgPb/genes/Pf_Pb_orth.xlsx')

xx.pf <- inner_join(xx.pf, pf.pb.orth, by = c('GeneID' = 'Pf'))
xx.pb <- inner_join(xx.pb, pf.pb.orth, by = c('GeneID' = 'Pb'))

venn.list <- list(Pf = xx.pf$Pb, Pb = xx.pb$GeneID)

library(ggVennDiagram)
g <- ggVennDiagram(venn.list, label_size = 6, set_size = 8,
                   set_color = c("Pb" = "firebrick","Pb" ="darkorchid3")) +  
  scale_color_manual(values = c("firebrick","darkorchid3")) +
  
  theme(legend.position = "none")


plot(g)


ggsave(filename="../Output/compScBdTgPb/figs/pf_pb_overlap.pdf", 
       plot=g,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)




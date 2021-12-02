library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
library(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(sctransform)
library(glassoFast)
library(igraph)
library(ggraph)
library(graphlayouts)
library(fdrtool)
library(parallel)

DNA.binding <- read.xlsx('../Input/compScBdTgPb/gene_function/DNA binding - modified from Data S2 Waldman-update 091321.xlsx', sheet = 1)
TGGT1.TGME49 <- read.xlsx('../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx')

DNA.binding <- left_join(DNA.binding, TGGT1.TGME49, by = c('ID' = 'TGME49'))

DNA.binding <- DNA.binding %>% dplyr::select(TGGT1, Description, Phenotype)
colnames(DNA.binding) <- c('GeneID', 'Description', 'Phenotype')

Kin <- read.xlsx('../Input/compScBdTgPb/gene_function/Tg GT1 kin&pho TXT search 090521.xlsx', sheet = 1)
Phos <- read.xlsx('../Input/compScBdTgPb/gene_function/Tg GT1 kin&pho TXT search 090521.xlsx', sheet = 2)

tg.prod.desc <- read.xlsx('../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx')
tg.opt.network <- readRDS('../Input/compScBdTgPb/RData/tg_opt_net.RData')
#tg.clust.net <- readRDS('../Input/compScBdTgPb/RData/tg_clust_net.RData')
#tg.clust.net.d <- readRDS('../Input/compScBdTgPb/RData/tg_clust_net_d.RData')
tg.opt.net.d <- readRDS('../Input/compScBdTgPb/RData/tg_opt_net_d.RData')
#tg.e.l <- readRDS('../Input/compScBdTgPb/RData/tg_e_l.RData')
tg.opt.e.l <- readRDS('../Input/compScBdTgPb/RData/tg_opt_e_l.RData')


DNA.binding <- left_join(DNA.binding, tg.opt.net.d, by = 'GeneID')
Kin <- left_join(Kin, tg.opt.net.d, by = c('Gene.ID' = 'GeneID'))
Phos <- left_join(Phos, tg.opt.net.d, by = c('Gene.ID' = 'GeneID'))


DNA.binding.present <- DNA.binding %>% dplyr::filter(!is.na(degree) )
DNA.binding.present <- DNA.binding.present %>% dplyr::select(GeneID, Description, degree, phase, Phenotype)
colnames(DNA.binding.present) <- c('GeneID', 'ProductDescription', 'degree', 'phase', 'Phenotype')

Kin.present <- Kin %>% dplyr::filter(!is.na(degree) )
Kin.present <- Kin.present %>% dplyr::select(Gene.ID, Product.Description, degree, phase)
Kin.present$Phenotype <- as.character('NA')
colnames(Kin.present) <- c('GeneID', 'ProductDescription', 'degree', 'phase', 'Phenotype')

Phos.present <- Phos %>% dplyr::filter(!is.na(degree))
Phos.present <- Phos.present %>% dplyr::select(Gene.ID, Product.Description, degree, phase)
Phos.present$Phenotype <- as.character('NA')
colnames(Phos.present) <- c('GeneID', 'ProductDescription', 'degree', 'phase', 'Phenotype')


DNA.binding.present$type <- 'DNA.binding'
Kin.present$type <- 'Kinase'
Phos.present$type <- 'Phosphatase'

DNA.binding.Kin.Phos.present <- rbind(DNA.binding.present, Kin.present, Phos.present) %>% arrange(desc(degree))

write.xlsx(DNA.binding.Kin.Phos.present, '../Output/compScBdTgPb/tables/tg_DNA_binding_Kin_Phos_network_degree_091321.xlsx') 


top.ROPs <- unique(DNA.binding.Kin.Phos.present$GeneID[grep('ROP', DNA.binding.Kin.Phos.present$ProductDescription)])[1:10]

ROP.net <- getNetFromVert(tg.clust.net, top.ROPs)
ROP.edge.list <- getEdgeList(ROP.net, tg.prod.desc)

write.xlsx(ROP7.edge.list, '../Output/compScBdTgPb/tables/tg_ROP7_edges.xlsx')

# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

p <- ggraph(ROP.net,layout = "stress") +
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = phase, size = 3),shape=21)+
  geom_node_text(aes(filter = size > 6, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme_graph(base_family = 'Helvetica') +
  theme(legend.position = "right")

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/bd_egress_clust4_network.pdf",
       plot=p4,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



# define a custom color palette
top.all <- unique(DNA.binding.Kin.Phos.present$GeneID)

all.net <- getNetFromVert(tg.clust.net, top.all)
all.edge.list <- getEdgeList(all.net, tg.prod.desc)


got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

p <- ggraph(all.net,layout = "stress") +
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = phase, size = 3),shape=21)+
  geom_node_text(aes(filter = size > 6, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme_graph(base_family = 'Helvetica') +
  theme(legend.position = "right")

plot(p)





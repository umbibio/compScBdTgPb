library(tidyverse)
library(bedtoolsr)
library(openxlsx)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(doParallel)
library(edgeR)



#######################################################
############# Peak Gene Assignment ####################
############ Raw counts per region ####################
#######################################################

gtf.file <- "../cut_and_run_pk/Genome/PlasmoDB-59_Pfalciparum3D7.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)

## using qval = 0.01 and gene size 69350000 bp in macs2 calling 
## narrow peaks called by macs2
narrow.peaks.file <-  "../cut_and_run_pk/Sample4-MK-Cut-20220903_S4.vs.Sample3-MK-Cut-20220903_S3_peaks.narrowPeak"
narrow.peaks <- read.table(narrow.peaks.file, header = F, sep = '\t', quote = NULL)


## on markov 
## cat *.narrowPeak | sort -k1,1 -k2,2n | mergeBed -i stdin > narrowPeaksUnion.bed


narrow.peaks$V4 <- paste(narrow.peaks$V1 , paste(narrow.peaks$V2, narrow.peaks$V3, sep= "-"), sep = ":")

peaks.all.sort <- narrow.peaks %>% arrange(V1, V2, V3)

## Get the transcripts
gtf.trans <- gtf %>% dplyr::filter(V3 == 'transcript')
parse.str <- strsplit(gtf.trans$V9, split = ' ')
inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
gtf.trans$V9 <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
gtf.trans$V9 <- gsub("\"", "", gtf.trans$V9)


## overlap the genes with peaks
genic.peaks <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = gtf.trans, wo = T)


## remaining peaks are intergenic
intergenic.peaks <- peaks.all.sort %>% dplyr::filter(!(V4 %in% genic.peaks$V4))
intergenic.peaks <- intergenic.peaks %>% arrange(V1, V2, V3) 

## Exons
gtf.exon <- gtf %>% dplyr::filter(V3 == 'exon')
gtf.exon.sort <- gtf.exon %>% arrange(V1, V4, V5)
parse.str <- strsplit(gtf.exon$V9, split = ' ')
inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
gtf.exon$V9 <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
gtf.exon$V9 <- gsub("\"", "", gtf.exon$V9)
gtf.exon <- gtf.exon %>% group_by(V9) %>% mutate(exon.ord = ifelse(V7 == '+', 1:n(), seq(n(), 1, by = -1)),
                                                 multiple.exon = ifelse(n() > 1, T, F))

## Filter for first exon coordinates (exon1 coordinates)
tmp.neg <- gtf.exon %>% filter(V7 == "-") %>% group_by(V9) %>%  dplyr::slice(which.max(V5))
tmp.pos <- gtf.exon %>% filter(V7 == "+") %>% group_by(V9) %>%  dplyr::slice(which.min(V5))
gtf.exon1 <- bind_rows(tmp.pos, tmp.neg)
gtf.exon1.sort <- gtf.exon1 %>% arrange(V1, V4, V5)

## Assign the peaks to nearest upstream gene (look at 5 closest in case of bi-directional)

intergenic.peaks.genes.dist <- bedtoolsr::bt.closest(a = intergenic.peaks, b = gtf.exon1.sort, D = "b", k = 5)

gene.id <- unlist(lapply(strsplit(intergenic.peaks.genes.dist$V19, split = ''), function(x){
  n = length(x)
  gene = paste(x[1:(n-1)], collapse = '')
  return(gene)
}))

intergenic.peaks.genes.dist$V19 <- gene.id 
## V16 <= 0 means the peak is at upstream 
## Find closest gene among top 5 that is upstreaam (min V16)
intergenic.peaks.genes.dist.upstream <- intergenic.peaks.genes.dist %>% dplyr::filter(V21 <= 0)
intergenic.peaks.genes.dist.upstream <- intergenic.peaks.genes.dist.upstream %>% group_by(V4) %>% 
  mutate(V22 = V21[which.min(abs(V21))])


## Filter the rest
intergenic.peaks.genes.dist.upstream <- intergenic.peaks.genes.dist.upstream %>% dplyr::filter(V21 == V22)

## filter the ones that are too far (1000 bp)
intergenic.peaks.genes.dist.upstream <- intergenic.peaks.genes.dist.upstream %>% dplyr::filter(abs(V21) < 1000)

intergenic.peaks <- intergenic.peaks.genes.dist.upstream %>% ungroup() %>%
  transmute(Chr = V1, strt = V2, stp = V3, peak_id = V4, GeneID = V19, type = 'itergenic')

gene.peaks <- genic.peaks %>% ungroup() %>%
  transmute(Chr = V1, strt = V2, stp = V3, peak_id = V4, GeneID = V19, type = 'within_gene')

all.peaks <- bind_rows(gene.peaks, intergenic.peaks)

write.xlsx(gene.peaks, "../cut_and_run_pk/all_peaks_gene_assignment.xlsx")



lactylated.genes.all <- all.peaks %>% group_by(GeneID) %>% transmute(peaks = list(peak_id), type = list(type)) %>% distinct()

lactylated.genes.within <- gene.peaks %>% group_by(GeneID) %>% transmute(peaks = list(peak_id), type = list(type)) %>% distinct()

lactylated.genes.intergenic <- intergenic.peaks %>% group_by(GeneID) %>% transmute(peaks = list(peak_id), type = list(type)) %>% distinct()


write.xlsx(lactylated.genes.all, "../cut_and_run_pk/lactylated_genes.xlsx")

# prep bed format 
# merge multiple peaks assigned to a single gene
# the duplicated peaks are the bidirectioonal peaks 

# peak.genes <- read.xlsx("../Input/Toxo_lab_adapt/BulkATACToxoPlasma/macs2_Union/Peaks_Genes_assigned.xlsx")
# peak.genes <- peak.genes %>% select(V1.x, V2.x, V3.x, V11, gene_name) 
# peak.genes.bed.merged <- peak.genes %>% arrange(V2.x) %>% 
#   group_by(gene_name) %>% mutate(start_peak = V2.x[which.min(V2.x)], end_peak = V3.x[which.max(V3.x)])  %>% 
#   mutate(V4 = ".", V5 = ".")
# 
# peak.genes.bed.merged.bed <- peak.genes.bed.merged %>% select(V1.x, start_peak, end_peak, V4, V5, V11, gene_name) %>%
#   distinct(gene_name, .keep_all = T)
# 
# 
# write.table(peak.genes.bed.merged.bed, "../Input/compScBdTgPb/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_merged_peaks_final.bed", 
#             sep = "\t", quote = F, row.names = F, col.names = F)
# write.xlsx(peak.genes.bed.merged.bed, "../Input/compScBdTgPb/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_merged_peaks_final.xlsx")

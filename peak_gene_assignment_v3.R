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

gtf.file <- "../Input/Toxo_genomics/Genome/ToxoDB-56_TgondiiGT1.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)

## using qval = 0.01 and gene size 69350000 bp in macs2 calling 
## narrow peaks called by macs2
narrow.peaks.dir <-  "../Input/Toxo_lab_adapt/BulkATACToxoPlasma/macs2_stringent/"
nr.peaks <- list.files(path = narrow.peaks.dir, pattern = ".narrowPeak")
nr.peaks <- nr.peaks[-8] # remove gDNA 

nr.peaks.file <- gsub("_", "-", nr.peaks)
x <- strsplit(nr.peaks.file,  "-")
samples <- c(sapply( x, "[", 5 )[-7], "RH")
indx <- sort(as.numeric(gsub("P", "", gsub("RH", "300",samples))),index.return =T)$ix
samples[indx]

nr.peaks.list <- list()
for(i in 1:length(nr.peaks)){
  
  in.dir <- paste(narrow.peaks.dir, nr.peaks[i], sep = "")
  in.file <- read.table(file = in.dir, header = F, sep = "", quote = "\"")
  #in.file <- in.file %>% filter(grepl('AAQM', V1))
  nr.peaks.list <- c( nr.peaks.list, list(in.file))
  
}
nr.peaks.list <- nr.peaks.list[indx]
names(nr.peaks.list) <- samples[indx]

## on markov 
## cat *.narrowPeak | sort -k1,1 -k2,2n | mergeBed -i stdin > narrowPeaksUnion.bed

narrowPeaksUnion <- read.table("../Input/Toxo_lab_adapt//BulkATACToxoPlasma/macs2_Union/narrowPeaksUnion.bed", 
                               header = F, sep = "\t", quote = NULL)


# generate narrow peaks for each sample using Union Peak regions file to load into IGV
out_dir <- "../Input/Toxo_lab_adapt/BulkATACToxo/Input/BulkATACToxoPlasma/macs2_Union"
UnionPeaksRegion.list <- lapply(1:length(nr.peaks.list), function(i){
  
  NAME = paste(paste(names(nr.peaks.list[i]), "Union", sep = "_"), "narrowPeak", sep = ".")
  sample.peaks <- nr.peaks.list[[i]]
  tmp <- bedtoolsr::bt.intersect(narrowPeaksUnion, sample.peaks, wo = T)
  tmp <- tmp %>% dplyr::select(!V4:V6) %>% group_by(V7) %>% dplyr::slice(which.max(V8))
  #write.table(tmp, paste(out_dir, NAME, sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)
  
  return(tmp)
})

names(UnionPeaksRegion.list) <- names(nr.peaks.list)


narrowPeaksUnion <- narrowPeaksUnion %>% filter(!grepl('AAQM*', V1))
narrowPeaksUnion$V4 <- paste(narrowPeaksUnion$V1 , paste(narrowPeaksUnion$V2, narrowPeaksUnion$V3, sep= "-"), sep = ":")

gtf.filt <- gtf %>% dplyr::filter(!grepl('AAQM*', V1))
peaks.all.sort <- narrowPeaksUnion %>% arrange(V1, V2, V3)



## Remove the first Exon from transcripts.
gtf.exon <- gtf.filt %>% dplyr::filter(V3 == 'exon')
gtf.exon.sort <- gtf.exon %>% arrange(V1, V4, V5)
parse.str <- strsplit(gtf.exon$V9, split = ' ')
inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
gtf.exon$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
gtf.exon$gene_name <- gsub("\"", "", gtf.exon$gene_name)
gtf.exon <- gtf.exon %>% group_by(V9) %>% mutate(exon.ord = ifelse(V7 == '+', 1:n(), seq(n(), 1, by = -1)),
                                                 multiple.exon = ifelse(n() > 1, T, F))
## Remove the exon1, but keep the Intron 1 , build exon2ton
gtf.exon.2Ton <- gtf.exon %>% mutate(V10 = ifelse(multiple.exon & V7 == '-', min(V4), 
                                                 ifelse(multiple.exon & V7 == '+', min(V5), V4)),
                                     V11 = ifelse(multiple.exon & V7 == '-', max(V4), 
                                                 ifelse(multiple.exon & V7 == '+', max(V5), V5))) %>%
  mutate(V4 = V10, V5 = V11) %>% 
  dplyr::select(-c(exon.ord,multiple.exon, V10, V11) ) %>% distinct()


## Overlap with peaks and filter peaks that are entirely within the genes.
## Overlapping peaks with exon2Ton and check to see if it is entirely within the gene. 
## Then from the sorted peaks we throw out all peaks entirely within the gene & not overlapping with exon1, 
## These peaks should not be assigned to any peaks.

peak.genes.ovrlp <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = gtf.exon.2Ton, wo = T)
peak.genes.filt <- peak.genes.ovrlp %>% dplyr::filter(V8  <= V2 & V9 >= V3)
peak.filt <- peaks.all.sort[!(peaks.all.sort$V4 %in%  peak.genes.filt$V4), ]
peak.filt.sort <- peak.filt %>% arrange(V1, V2, V3)


## filter gtf for transcripts only to get the coordinates of start and end of gene
gtf.filt.trn <- gtf.filt %>% filter(V3 == "transcript")
gtf.filt.trn$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt.trn$V9)))
gtf.filt.trn$gene_name <- gsub("\"", "", gtf.filt.trn$gene_name)

## Filter for first exon coordinates (exon1 coordinates)
tmp.neg <- gtf.exon %>% filter(V7 == "-") %>% group_by(V9) %>%  dplyr::slice(which.max(V5))
tmp.pos <- gtf.exon %>% filter(V7 == "+") %>% group_by(V9) %>%  dplyr::slice(which.min(V5))
gtf.exon1 <- bind_rows(tmp.pos, tmp.neg)
gtf.exon1.sort <- gtf.exon1 %>% arrange(V1, V4, V5)

## Assign the peaks to nearest upstream gene (look at 5 closest in case of bi-directional)

peaks.genes.dist <- bedtoolsr::bt.closest(a = peak.filt.sort, b = gtf.exon1.sort, D = "b", k = 5)
peaks.genes.dist$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", peaks.genes.dist$V13)))

peaks.genes.dist.trns <- left_join(peaks.genes.dist, gtf.filt.trn, by = "gene_name")

## V16 is the distance of the peak to the exon 1 
## we need to overcome the issue with the  ones with  dist = 0
## on pos strand V3.x (end of peak) should not exceed V5.y (end of transcript/exon_n)
## on neg strand V2.x (start of peak) is not less than V4.y (beggining of the transcript/exon_1)

peaks.genes.dist.trns <- peaks.genes.dist.trns %>% filter(!(V16 == 0 & V11 == "+" & V3.x > V5.y))

peaks.genes.dist.trns <- peaks.genes.dist.trns %>% filter(!(V16 == 0 & V11 == "-" & V2.x < V4.y))


## V16 <= 0 means the peak is at upstream 
## Find closest gene among top 5 that is upstreaam (min V16)
peaks.genes.dist.trns <- peaks.genes.dist.trns %>% dplyr::filter(V16 <= 0)
peaks.genes.dist.trns <- peaks.genes.dist.trns %>% group_by(V4.x) %>% 
  mutate(V17 = V16[which.min(abs(V16))])


## Filter the rest
peaks.genes.dist.trns <- peaks.genes.dist.trns %>% dplyr::filter(V16 == V17)
# parse.str <- strsplit(peaks.genes.dist$V13, split = ';')
# inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
# peaks.genes.dist$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))


## filter the ones that are too far (3000 bp)
peaks.genes.dist.trns.filt <- peaks.genes.dist.trns %>% dplyr::filter(abs(V16) < 3000)

tmp <- peaks.genes.dist.trns.filt %>% group_by(V4.x) %>% summarise(total = length(unique(gene_name)))
tmp.2 <- peaks.genes.dist.trns.filt %>% group_by(gene_name) %>% summarise(total = length(unique(V4.x)))

write.xlsx(peaks.genes.dist.trns.filt, "../Input/Toxo_lab_adapt/BulkATACToxoPlasma/macs2_Union/Peaks_Genes_assigned.xlsx")

# prep bed format 
# merge multiple peaks assigned to a single gene
# the duplicated peaks are the bidirectioonal peaks 

peak.genes <- read.xlsx("../Input/Toxo_lab_adapt/BulkATACToxoPlasma/macs2_Union/Peaks_Genes_assigned.xlsx")
peak.genes <- peak.genes %>% select(V1.x, V2.x, V3.x, V11, gene_name) 
peak.genes.bed.merged <- peak.genes %>% arrange(V2.x) %>% 
  group_by(gene_name) %>% mutate(start_peak = V2.x[which.min(V2.x)], end_peak = V3.x[which.max(V3.x)])  %>% 
  mutate(V4 = ".", V5 = ".")

peak.genes.bed.merged.bed <- peak.genes.bed.merged %>% select(V1.x, start_peak, end_peak, V4, V5, V11, gene_name) %>%
  distinct(gene_name, .keep_all = T)


write.table(peak.genes.bed.merged.bed, "../Input/compScBdTgPb/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_merged_peaks_final.bed", 
            sep = "\t", quote = F, row.names = F, col.names = F)
write.xlsx(peak.genes.bed.merged.bed, "../Input/compScBdTgPb/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_merged_peaks_final.xlsx")



##########################



## generate counts for each sample using Union peak regions from sorted bam files of each sample
## convert peakRegion.bed to saf format acceptable by feature count 
## awk 'OFS="\t" {print $1"_"$2"_"$3, $1, $2, $3, "."}' peaksRegion.bed > peaksRegion.saf
##  screen
## bsub -Is -q interactive -n5 -R "span[hosts=1]" /bin/bash
## module load subread/1.6.2
## featureCounts  -a peaksRegion.saf -F SAF -p -o peaks_count.txt /project/umb_kourosh_zarringhalam/yasaman/ATACBulkToxo/sam/*.sorted.bam 




## add gene names to peak region in the count table

count <- read.table("../Input/Toxo_lab_adapt/BulkATACToxoPlasma/macs2_Union/peaks_count.txt", header = T, sep = '\t', stringsAsFactors = F)
count$Geneid <- sub("^([^_]+)_([^_]+)_([^_]+)_", "\\1_\\2:\\3-", count$Geneid)
count <- count %>% dplyr::select(Geneid, colnames(count)[grepl("P[0-9]|RH", colnames(count))])


Peaks_Genes_assigned <- read.xlsx("../Input/Toxo_lab_adapt/BulkATACToxoPlasma/macs2_Union/Peaks_Genes_assigned.xlsx")
Peaks_Genes_assigned <- Peaks_Genes_assigned %>% dplyr::transmute(gene_name =  gene_name, peak_location = V4.x)

peak.gene.count <- left_join(Peaks_Genes_assigned, count, by = c( "peak_location" = "Geneid"))  
peak.gene.count.ave <- peak.gene.count %>% select(-peak_location) %>% group_by(gene_name) %>% summarise_all("mean")
peak.gene.count.ave <- peak.gene.count.ave %>% group_by(gene_name) %>% summarise_all("ceiling")


# rename samples (column)
X <- lapply(strsplit(colnames(peak.gene.count.ave), "\\_"), `[[`, 1)
XX <- sub(".*\\.", "", X) 
XX <- gsub("gene", "gene_name", XX)

colnames(peak.gene.count.ave) <- XX

write.xlsx(peak.gene.count.ave, "../Input/Toxo_lab_adapt/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_raw_counts_avg_final.xlsx")



## up to here ##

# narrow.peaks.dir <-  "../Input/Toxo_lab_adapt/BulkATACToxoPlasma/macs2_stringent/"
# nr.peaks <- list.files(path = narrow.peaks.dir, pattern = ".narrowPeak")
# 
# 
# 
# Names <- gsub("_peaks", "", sapply(X, "[", 1))
# cond <- rep("extra", 8) 
# passage <- gsub("_S[0-9]","", c(sapply(strsplit(Names, "-"), "[" , 5)[-c(7,8)], "RH","RH_gDNA"))
# treatment <- paste(cond, passage, sep = ".")
# Sample <- gsub(".*_", "", Names)
# 
# expmnt <- data.frame(Names, cond, passage, treatment, Sample)
# expmnt
# 
# colnames(peak.gene.count)[3:10] <- treatment
# 
# write.xlsx(peak.gene.count, "../Input/compScBdTgPb/BulkATACToxo/Input/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_raw_counts.xlsx")
# 

# #######################################################
# ################ Normaliization #######################
# #######################################################
# 
# count.tab <- read.xlsx("../Input/compScBdTgPb/BulkATACToxo/Input/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_raw_counts.xlsx")
# 
# # cpm normalization
# CPM <- cbind(count.tab[,1:2],cpm(count.tab[,3:ncol(count.tab)]))
# CPM.avg <- CPM[,-2] %>% group_by(gene_name) %>% summarise_all("mean")
# 
# # filter low expressed genes
# keep <- rowSums(CPM.avg > 2) >= 3
# 
# 
# count.avg <- count[,-2] %>% group_by(gene_name) %>% summarise_all("mean")
# count.avg <-cbind(gene_id = count.avg$gene_id, round(count.avg[,-1], digits = 0))
# #count.avg <- count.avg %>% dplyr::select(-contains("gDNA"))
# 
# x <- count.avg[keep, 2:ncol(count.avg)]
# rownames(x) <- count.avg$gene_name[keep]
# treatment <- expmnt$treatment
# treatment <- factor(treatment, levels = unique(treatment))
# 
# 
# # Creating DGElist object
# y <- DGEList(counts = x, group = treatment) 
# y <- calcNormFactors(y)
# y$samples
# 
# # logarithmic scale
# logCPM <- cpm(y, log=TRUE, prior.count=3, normalized.lib.sizes=TRUE) 
# plotMDS(logCPM, cex = 0.8)
# 
# logCPM.df <- logCPM %>% data.frame() %>% rownames_to_column(var = "GeneID")
# logCPM.df <- logCPM.df %>% mutate(across(c(2:8),.fns = ~./extra.RH_gDNA))
# logCPM.df <- logCPM.df %>% dplyr::select(!contains("gDNA"))
# write.xlsx(logCPM.df, "../Input/BulkATACToxoPlasma/macs2_Union/Toxo_tab_ATAC_seq_logCPM.xlsx")
# 
# 

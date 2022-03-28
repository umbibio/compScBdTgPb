library(tidyverse)
library(bedtoolsr)
#library(DiffBind)
library(openxlsx)

gtf.file <- "../Input/compScBdTgPb/BulkATACToxo/Input/BulkATACToxoPlasma/Genome/ToxoDB-56_TgondiiGT1.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)


# finding first exons according to the strand of gene ## No need to run this

# gtf.df <- gtf %>% dplyr::select(c(V3:V5, V7,V9))
# str <- strsplit(gtf.df$V9, " ")
# gtf.df$V9 <- gsub(";", "",gsub("\"", "", lapply(str, "[" , 4)))
# gtf.df <- gtf.df %>% dplyr::select(V9, everything())
# gtf.exon <- gtf.df %>% filter(V3 == "exon") 
# tmp.neg <- gtf.exon %>% group_by(V9) %>% filter(V7 == "-") %>% dplyr::slice(which.max(V5))
# tmp.pos <- gtf.exon %>% group_by(V9) %>% filter(V7 == "+") %>% dplyr::slice(which.min(V4))
# 
# exon <- rbind (tmp.neg, tmp.pos)
# colnames(exon) <- c("GeneName", "first_ex", "first_ex_start", "first_ex_stop", "strand")
# 
# gtf.transcripit <- gtf.df %>% filter(V3 == "transcript")
# colnames(gtf.transcripit) <- c("GeneName", "Name", "gene_start", "gene_stop", "strand")
# 
# trans.first.exon <- left_join(gtf.transcripit, exon, by = c("GeneName", "strand"))


## using qval = 0.01 and gene size 69350000 bp in macs2 calling 
## narrow peaks called by macs2
narrow.peaks.dir <-  "../Input/compScBdTgPb/BulkATACToxo/Input/BulkATACToxoPlasma/macs2_stringent/"
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

narrowPeaksUnion <- read.table("../Input/compScBdTgPb/BulkATACToxo/Input/BulkATACToxoPlasma/macs2_Union/narrowPeaksUnion.bed", 
                               header = F, sep = "\t", quote = NULL)


# generate narrow peaks for each sample using Union Peak regioins file to load into IGV
out_dir <- "../Input/compScBdTgPb/BulkATACToxoPlasma/macs2_Union"
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


##

## Remove the first Exon from transcripts.
gtf.exon <- gtf.filt %>% dplyr::filter(V3 == 'exon')
gtf.exon.sort <- gtf.exon %>% arrange(V1, V4, V5)
parse.str <- strsplit(gtf.exon$V9, split = ' ')
inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
gtf.exon$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
gtf.exon$gene_name <- gsub("\"", "", gtf.exon$gene_name)
gtf.exon <- gtf.exon %>% group_by(V9) %>% mutate(exon.ord = ifelse(V7 == '+', 1:n(), seq(n(), 1, by = -1)),
                                                 multiple.exon = ifelse(n() > 1, T, F))
## Remove the exon1, but keep the Intron 1
gtf.exon.2Ton <- gtf.exon %>% mutate(V10 = ifelse(multiple.exon & V7 == '-', min(V4), 
                                                 ifelse(multiple.exon & V7 == '+', min(V5), V4)),
                                     V11 = ifelse(multiple.exon & V7 == '-', max(V4), 
                                                 ifelse(multiple.exon & V7 == '+', max(V5), V5))) %>%
  mutate(V4 = V10, V5 = V11) %>% 
  dplyr::select(-c(exon.ord,multiple.exon, V10, V11) ) %>% distinct()


## Overlap with peaks and filter peaks that are entirely within the genes.
peak.genes.ovrlp <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = gtf.exon.2Ton, wo = T)
peak.genes.filt <- peak.genes.ovrlp %>% dplyr::filter(V8  <= V2 & V9 >= V3)
peak.filt <- peaks.all.sort[!(peaks.all.sort$V4 %in%  peak.genes.filt$V4), ]
peak.filt.sort <- peak.filt %>% arrange(V1, V2, V3)


## Filter for first exon coordinates
tmp.neg <- gtf.exon %>% filter(V7 == "-") %>% group_by(V9) %>%  dplyr::slice(which.max(V5))
tmp.pos <- gtf.exon %>% filter(V7 == "+") %>% group_by(V9) %>%  dplyr::slice(which.min(V5))
gtf.exon1 <- bind_rows(tmp.pos, tmp.neg)
gtf.exon1.sort <- gtf.exon1 %>% arrange(V1, V4, V5)

## Assign the peaks to nearest upstream gene (look at 5 closest in case of bi-directional)
peaks.genes.dist <- bedtoolsr::bt.closest(a = peak.filt.sort, b = gtf.exon1.sort, D = "b", k = 5)

## Find closest gene among top 5 that is upstreaam
peaks.genes.dist <- peaks.genes.dist %>% dplyr::filter(V16 <= 0)
peaks.genes.dist <- peaks.genes.dist %>% group_by(V4) %>% 
  mutate(V17 = V16[which.min(abs(V16))])


## Filter the rest
peaks.genes.dist <- peaks.genes.dist %>% dplyr::filter(V16 == V17)
parse.str <- strsplit(peaks.genes.dist$V13, split = ';')
inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
peaks.genes.dist$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))




############### STOPPED HERE

peaks.genes.dist <- peaks.genes.dist %>% ungroup() %>% transmute(Chr = V1, gene_start = V8, gene_stop = V9,
                                                                 gene_lenght = V9 - V8 ,
                                                                 strand = V11, gene_id = gene_name, peak_location = V4, 
                                                                 peak_start = V2, peak_stop = V3, distance = V14)

peaks.intergenic.sort.filt <- peaks.intergenic.sort[peaks.intergenic.sort$V4 %in% peaks.genes.dist$peak_location,]
# Find the ones that overlap and calculate the percantage off overlap
peak.genes.ovrlp <- bedtoolsr::bt.intersect(a = peaks.intergenic.sort.filt, b = gtf.trans.sort, wo = T)
parse.str <- strsplit(peak.genes.ovrlp$V13, split = ' ')
inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
peak.genes.ovrlp$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
peak.genes.ovrlp$gene_name <- gsub("\"", "", peak.genes.ovrlp$gene_name)

peak.genes.ovrlp <- peak.genes.ovrlp %>% ungroup() %>% transmute(gene_id = gene_name,
                                                                 peak_location = V4, peak_start = V2, 
                                                                 peak_stop = V3, overlap = V14, 
                                                                 overlap_pct = V14 / (V9 - V8))
peak.genes.dist.ovlp <- left_join(peaks.genes.dist,peak.genes.ovrlp, 
                                  by = c( "gene_id","peak_location", "peak_start" , "peak_stop"))
peak.genes.dist.ovlp$overlap_pct[is.na(peak.genes.dist.ovlp$overlap_pct)] <- 0
peak.genes.dist.ovlp$overlap[is.na(peak.genes.dist.ovlp$overlap)] <- 0

# to add the flag for upstream and downstream if needed
peak.genes.dist.ovlp <- peak.genes.dist.ovlp %>% mutate(Peak_gene_start_dist = peak_start - gene_start)
peak.genes.dist.ovlp <- peak.genes.dist.ovlp %>% mutate(peak_gene_end_dist = peak_stop - gene_stop)


peak.genes.dist.ovlp <- peak.genes.dist.ovlp %>% mutate(direction = case_when(
  (strand == "-" & peak_gene_end_dist > 0) ~ "upstream",
  (strand == "-" & peak_gene_end_dist < 0) ~ "downstream",
  (strand == "+" & Peak_gene_start_dist < 0) ~ "upstream",
  (strand == "+" & Peak_gene_start_dist > 0) ~ "downstream"))



## filtering 
peaks.genes.filt <- peak.genes.dist.ovlp %>% 
  dplyr::filter(((strand == '+' & distance >= 0) | (strand == '-' & distance <= 0)) & 
                  (abs(distance) < 3000)  & overlap_pct <= 0.5 
                #direction == "upstream"
  )
write.xlsx(peaks.genes.filt, "../Input/BulkATACToxoPlasma/macs2_Union/Peaks_Genes_filtered.xlsx")

peaks.Region <- peaks.genes.filt %>% 
  dplyr::select(Chr,  peak_start, peak_stop, peak_location) %>% distinct(peak_location,  .keep_all = T)
write.table(peaks.Region, file="../Input/BulkATACToxoPlasma/macs2_Union/peaksRegion.bed",
            quote=F, sep="\t", row.names=F, col.names=F)


## generate counts for each sample using Union peak regions from sorted bam files of each sample
## convert peakRegion.bed to saf format acceptable by feature count 
## awk 'OFS="\t" {print $1"_"$2"_"$3, $1, $2, $3, "."}' peaksRegion.bed > peaksRegion.saf
##  screen
## bsub -Is -q interactive -n5 -R "span[hosts=1]" /bin/bash
## module load subread/1.6.2
## featureCounts  -a peaksRegion.saf -F SAF -p -o peaks_count.txt /project/umb_kourosh_zarringhalam/yasaman/ATACBulkToxo/sam/*.sorted.bam 




## add gene names to peak region count table

# peaks.genes.filt <- read.xlsx("../Input/BulkATACToxoPlasma/macs2_Union/Peaks_Genes_filtered.xlsx")

count <- read.table("../Input/BulkATACToxoPlasma/macs2_Union/peaks_count.txt", header = T, sep = '\t', stringsAsFactors = F)
count$Geneid <- sub("^([^_]+)_([^_]+)_([^_]+)_", "\\1_\\2:\\3-", count$Geneid)
count <- count %>% dplyr::select(Geneid, colnames(count)[grepl("P[0-9]|RH", colnames(count))])


peak.gene.count <- left_join(peaks.genes.filt, count, by = c( "peak_location" = "Geneid"))  
peak.gene.count <- peak.gene.count %>% dplyr::select(gene_id, peak_location, colnames(count)[-1]) 


narrow.peaks.dir <-  "../Input/BulkATACToxoPlasma/macs2_stringent/"
nr.peaks <- list.files(path = narrow.peaks.dir, pattern = ".narrowPeak")
X <- strsplit(nr.peaks, "\\.")

Names <- gsub("_peaks", "", sapply(X, "[", 1))
cond <- rep("extra", 8) 
passage <- gsub("_S[0-9]","", c(sapply(strsplit(Names, "-"), "[" , 5)[-c(7,8)], "RH","RH_gDNA"))
treatment <- paste(cond, passage, sep = ".")
Sample <- gsub(".*_", "", Names)

expmnt <- data.frame(Names, cond, passage, treatment, Sample)
expmnt

colnames(peak.gene.count)[3:10] <- treatment

write.xlsx(peak.gene.count, "../Input/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_raw_counts.xlsx")







# this has not been done (merge)

## merge multiple peaks assigned to a sinngle gene 
#peaks.genes.filt.merged <- peaks.genes.filt %>% group_by(gene_id) %>% 
# mutate(min_peak_start = min(peak_start), max_peak_stop = max(peak_stop))
# peaks.genes.filt.merged <- peaks.genes.filt %>% group_by(Chr, gene_id) %>% 
#   summarise(min_peak_start = min(peak_start), max_peak_stop = max(peak_stop)) 
# 
# peaks.genes.filt.merged$merged_peak_location <- paste(peaks.genes.filt.merged$Chr, 
#                                                paste(peaks.genes.filt.merged$min_peak_start, 
#                                                      peaks.genes.filt.merged$max_peak_stop, sep = "-"), sep =":")
# 
# peaks.Region.merged <- peaks.genes.filt.merged %>% 
#   dplyr::select(Chr,  min_peak_start, max_peak_stop, merged_peak_location)
# colnames(peaks.Region.merged) <- c('chrom', 'chromStart', 'chromEnd')
# 
# 
# write.table(peaks.Region, file="../Input/BulkATACToxoPlasma/macs2_Union/peaksRegionMerged.bed",
#             quote=F, sep="\t", row.names=F, col.names=F)
# 



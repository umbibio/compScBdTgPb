library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(doParallel)


source('./util_funcs.R')

## IDs
prod.desc  <- read.xlsx('../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx')



## Read in the SNP file from MJ's 2013 BMC Genomic publication (take the first 12 columns)
SNP.syn <- read.xlsx('../Input/compScBdTgPb/RH_SNPs/SNPs_BMC_Genomics_KZ.xlsx', sheet = 1)
SNP.nonsyn <- read.xlsx('../Input/compScBdTgPb/RH_SNPs/SNPs_BMC_Genomics_KZ.xlsx', sheet = 2)
INDEL.CDS <- read.xlsx('../Input/compScBdTgPb/RH_SNPs/SNPs_BMC_Genomics_KZ.xlsx', sheet = 3)
SNP.noCDS <- read.xlsx('../Input/compScBdTgPb/RH_SNPs/SNPs_BMC_Genomics_KZ.xlsx', sheet = 4)
INDEL.noCDS <- read.xlsx('../Input/compScBdTgPb/RH_SNPs/SNPs_BMC_Genomics_KZ.xlsx', sheet = 5)
SNP.intergenic <- read.xlsx('../Input/compScBdTgPb/RH_SNPs/SNPs_BMC_Genomics_KZ.xlsx', sheet = 6)
INDEL.intergenic <- read.xlsx('../Input/compScBdTgPb/RH_SNPs/SNPs_BMC_Genomics_KZ.xlsx', sheet = 7)

my.col.names <- c("Chromosome", "Position", "Depth", "GT1_ref", "RH_SNP", 
                  "Flanking", "GeneID", "ProductDescription", "CodonChange", "AAChange", 
                  "Strand", "GeneType")

## take the first 12 columns
SNP.syn <- SNP.syn[,1:12]
colnames(SNP.syn) <- my.col.names
SNP.syn$SNPIndel <- "SNP"
SNP.syn$Type <- "syn"
nrow(SNP.syn)

SNP.nonsyn <- SNP.nonsyn[,1:12]
colnames(SNP.nonsyn) <- my.col.names
SNP.nonsyn$SNPIndel <- "SNP"
SNP.nonsyn$Type <- "nonsyn"
nrow(SNP.nonsyn)

INDEL.CDS <- INDEL.CDS[,1:12]
colnames(INDEL.CDS) <- my.col.names
INDEL.CDS$SNPIndel <- "INDEL"
INDEL.CDS$Type <- "cds"
nrow(INDEL.CDS)

SNP.noCDS <- SNP.noCDS[,1:12]
colnames(SNP.noCDS) <- my.col.names
SNP.noCDS$SNPIndel <- "SNP"
SNP.noCDS$Type <- "nocds"
nrow(SNP.noCDS)

INDEL.noCDS <- INDEL.noCDS[,1:12]
colnames(INDEL.noCDS) <- my.col.names
INDEL.noCDS$SNPIndel <- "INDEL"
INDEL.noCDS$Type <- "nocds"
nrow(INDEL.noCDS)

SNP.intergenic <- SNP.intergenic[,1:12]
colnames(SNP.intergenic) <- my.col.names
SNP.intergenic$SNPIndel <- "SNP"
SNP.intergenic$Type <- "intergenic"
nrow(SNP.intergenic)

INDEL.intergenic <- INDEL.intergenic[,1:12]
colnames(INDEL.intergenic) <- my.col.names
INDEL.intergenic$SNPIndel <- "INDEL"
INDEL.intergenic$Type <- "intergenic"

SNPs <- rbind(SNP.nonsyn, SNP.syn, INDEL.CDS, SNP.noCDS, INDEL.noCDS, SNP.intergenic, INDEL.intergenic)
SNPs.stats <- SNPs %>% group_by(SNPIndel, Type) %>% summarise(total = n())

p <-  ggplot(SNPs.stats, aes(x=SNPIndel, y=total, fill = Type)) +
  geom_bar(stat="identity", position = position_dodge(width = .9))+
  geom_text(aes(label=total, group = Type), position = position_dodge(width = .9),
            vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  ylab('Total') + xlab('') +
  #scale_x_discrete(labels=c("Intra", "Extra", "Brady")) +
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


ggsave(filename="../Output/compScBdTgPb/figs/SNP_totals.pdf", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## First, check if any Bradyzoite/Tachyzoite genes have SNP/Indel associated with them
Intra.vs.Brady.markers.sig <- read.xlsx('../Output/compScBdTgPb/tables/Table3_Intra_vs_brady_markers.xlsx')
Brady.genes <- Intra.vs.Brady.markers.sig[Intra.vs.Brady.markers.sig$avg_log2FC < 0, ]
any(Brady.genes$GeneID %in% SNPs$GeneID) ## Does not look like it.
Intra.genes <- Intra.vs.Brady.markers.sig[Intra.vs.Brady.markers.sig$avg_log2FC > 0, ]
any(Intra.genes$GeneID %in% SNPs$GeneID) ## Does not look like it.

## Wrtie a fasta file
Xfasta <- character(nrow(SNPs) * 2)
Xfasta[c(TRUE, FALSE)] <- paste(paste0(">", SNPs$Chromosome), 
                                SNPs$Position, SNPs$GT1_ref, SNPs$RH_SNP, SNPs$GeneID, sep = '_')
Xfasta[c(FALSE, TRUE)] <- SNPs$Flanking
writeLines(Xfasta, "../Input/compScBdTgPb/RH_SNPs/flanking.fasta")


### Aling the above reads to Toxo Genome V 56
# bwa index ToxoDB-56_TgondiiGT1_Genome.fasta
# bwa mem -t 8 ../genome/ToxoDB-56_TgondiiGT1_Genome.fasta flanking.fasta > flanking.sam
# bamToBed -i flanking.sam > flanking.bed 
# bedtools getfasta -s -fi ../genome/ToxoDB-56_TgondiiGT1_Genome.fasta -bed flanking.bed -fo flanking_GT1.fasta
# Sequences in flanking_GT1.fasta and flanking.fasta match exactly. The predicted SNP at the middle seems to be accurate
# cat ToxoDB-56_TgondiiGT1.gtf | awk '{if ($3 == "transcript") print $0}' > ToxoDB-56_TgondiiGT1_transcripts.gtf
# bedtools sort -i flanking.bed > flanking_sorted.bed
# bedtools sort -i ToxoDB-56_TgondiiGT1_transcripts.gtf ToxoDB-56_TgondiiGT1_transcripts_sorted.gtf
# bedtools closest -D "b" -a flanking_sorted.bed -b ../genome/ToxoDB-56_TgondiiGT1_transcripts_sorted.gtf > flanking_closest.bed

flanking_closest <- read.table('../Input/compScBdTgPb/RH_SNPs/flanking_closest.bed',
                               sep = '\t', header = F, quote = '#')


# Identify flanking reagions with BFD1 binding site
BFD1.site <- "(A|G)CTGG"
BFD1.site.rc <- "CCAG(C|T)"

potential.match <- gregexec(BFD1.site, SNPs$Flanking)
potential.match.rc <- gregexec(BFD1.site.rc, SNPs$Flanking)

positions <- unlist(lapply(potential.match, `[[`, 1))
SNPs$BFD1.site <- positions
SNPs$BFD1.site.overlaps <- ifelse(positions > 43 & positions < 51, 'yes', 'no')

positions.rc <- unlist(lapply(potential.match.rc, `[[`, 1))
SNPs$BFD1.site.rc <- positions.rc
SNPs$BFD1.site.rc.overlaps <- ifelse(positions.rc > 43 & positions.rc < 51, 'yes', 'no')

SNPs$ID <- paste(SNPs$Chromosome, SNPs$Position, SNPs$GT1_ref, SNPs$RH_SNP, SNPs$GeneID, sep = '_')

SNPs <- left_join(SNPs, flanking_closest, by = c('ID' = 'V4'))


## Read in the BFD1 binding sites
BFD1.cut_and_run <- read.xlsx('../Input/compScBdTgPb/RH_SNPs/BFD1_cut_and_run_KZ.xlsx')
BFD1.cut_and_run$Gene1.ID <- TGGT1_ME49$TGGT1[match(BFD1.cut_and_run$Gene1_ID, TGGT1_ME49$TGME49)]
BFD1.cut_and_run$Gene2.ID <- TGGT1_ME49$TGGT1[match(BFD1.cut_and_run$Gene2_ID, TGGT1_ME49$TGME49)]
SNPs$V17 <- gsub(" ", "", gsub("\"", "", gsub('gene_id', '', gsub('.*;', '', SNPs$V15))))

SNPs$cut_run_gene1 <- BFD1.cut_and_run$Gene1.ID[match(SNPs$V17, BFD1.cut_and_run$Gene1.ID)] 
SNPs$cut_run_gene1_motif <- BFD1.cut_and_run$Gene1_has_motif[match(SNPs$V17, BFD1.cut_and_run$Gene1.ID)] 
SNPs$cut_run_gene1_reg <- BFD1.cut_and_run$Gene1_group[match(SNPs$V17, BFD1.cut_and_run$Gene1.ID)] 
SNPs$cut_run_gene1_peak_motif <- BFD1.cut_and_run$Peak_has_motif[match(SNPs$V17, BFD1.cut_and_run$Gene1.ID)] 

SNPs$cut_run_gene2 <- BFD1.cut_and_run$Gene2.ID[match(SNPs$V17, BFD1.cut_and_run$Gene2.ID)] 
SNPs$cut_run_gene2_motif <- BFD1.cut_and_run$Gene2_has_motif[match(SNPs$V17, BFD1.cut_and_run$Gene2.ID)] 
SNPs$cut_run_gene2_reg <- BFD1.cut_and_run$Gene2_group[match(SNPs$V17, BFD1.cut_and_run$Gene2.ID)] 
SNPs$cut_run_gene2_peak_motif <- BFD1.cut_and_run$Peak_has_motif[match(SNPs$V17, BFD1.cut_and_run$Gene2.ID)] 

SNPs <- left_join(SNPs, prod.desc, by = c('V17' = 'GeneID'))
write.xlsx(SNPs, '../Input/compScBdTgPb/RH_SNPs/SNP_BFD1_overlap.xlsx')


SNPs.filt <- SNPs %>% dplyr::filter((BFD1.site != -1 & !is.na(cut_run_gene1)) | 
                                      (BFD1.site != -1 & !is.na(cut_run_gene2)) | 
                                      (BFD1.site.rc != -1 & !is.na(cut_run_gene1)) |
                                      (BFD1.site.rc != -1 & !is.na(cut_run_gene2)))

SNPs.candides <- SNPs.filt %>% transmute(GT1_ref = GT1_ref, RH_SNP = RH_SNP, CodonChange = CodonChange,
                                         AAChange = AAChange, Flanking = Flanking, 
                                         BFD1.site = BFD1.site, BFD1.site.rc = BFD1.site.rc, 
                                         nearest_gene = V17, distnace_to_nearest_gene = V16,
                                         ProductDescription = ProductDescription.y,
                                         cut_run_gene1 = cut_run_gene1_reg,
                                         cut_run_gene2 = cut_run_gene2_reg)

write.xlsx(SNPs.candides, '../Input/compScBdTgPb/RH_SNPs/SNP_BFD1_candidates.xlsx')

gene.id <- 'TGGT1_312100'
gene.id %in% Brady.genes$GeneID
gene.id %in% Intra.genes$GeneID



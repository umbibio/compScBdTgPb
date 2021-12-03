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
library(glassoFast)
library(igraph)
library(ggraph)
library(graphlayouts)
library(Signac)
library(Seurat)
library(patchwork)
library(hdf5r)
library(GenomeInfoDbData)
library(GenomicRanges)
library(GenomicAlignments)
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)

ME49.fasta <- readDNAStringSet("../Input/ME49/ToxoDB-52_TgondiiME49_Genome.fasta")
chrs <- names(ME49.fasta)[grep("TGME49_chr", names(ME49.fasta))]

chr.len <- data.frame(chr = gsub(" ", "", unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 1))),
                      len = as.numeric(gsub('length=', '', unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 4)))))

txdb <- makeTxDbFromGFF(file="../Input/ME49/ToxoDB-52_TgondiiME49_filter.gtf",
                        dataSource="Toxodb",
                        organism="Toxoplasma")

trans_biotypes <- select(txdb, keys=keys(txdb, "TXID"), 
                         columns = "TXTYPE", keytype =  "TXID")


genome(txdb) <- 'ME49'
tx_genes <- genes(txdb)
tmp <- chr.len$len[match(names(seqlengths(txdb)), chr.len$chr)]
names(tmp) <- names(seqlengths(txdb))
seqlengths(tx_genes) <- tmp

seqlevels(tx_genes)
inds <- c(5,6,1,2,3,7,8,10,11,9,4,12,13,14)

#seqlevels(tx_genes) <- gsub('TGME49_', '', seqlevels(tx_genes)[inds])
seqlevels(tx_genes) <- seqlevels(tx_genes)[inds]
isCircular(tx_genes) <- rep(F, length(isCircular(tx_genes)))

seqinfo(tx_genes)

counts <- Read10X_h5(filename = "../Input/scATAC/ME49_cell_ranger/raw_peak_bc_matrix.h5")
peak_anno <- read_tsv("../Input/scATAC/ME49_cell_ranger/raw_peak_bc_matrix/peaks.bed")

counts <- Read10X_h5(filename = "../Input/scATAC/ME49_cell_ranger/filtered_peak_bc_matrix.h5")
peak_anno <- read_tsv("../Input/scATAC/ME49_cell_ranger/peak_annotation.tsv")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = seqinfo(tx_genes),
  fragments = '../Input/scATAC/ME49_cell_ranger/fragments.tsv.gz',
  min.cells = 3,
  min.features = 10
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)


pbmc[['peaks']]
granges(pbmc)
annotations <- tx_genes

#seqlevelsStyle(annotations) <- 'UCSC'

annotations$gene_biotype <- 'protein_coding'
annotations$gene_name <- annotations$gene_id
Annotation(pbmc) <- annotations
pbmc <- NucleosomeSignal(object = pbmc)
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

peaks <- CallPeaks(
  object = pbmc,
  macs2.path = "/Users/kouroshz/miniconda3/envs/macs2/bin/macs2"
)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc


pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

DepthCor(pbmc)


pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(pbmc)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c("TGME49-292950","TGME49-293180","TGME49-293258"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

peaks <- CallPeaks(
  object = pbmc,
  group.by = "predicted.id",
  macs2.path = "/home/stuartt/miniconda3/envs/signac/bin/macs2"
)



require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)


library(rtracklayer)
library(GenomicRanges)




# peak-barcode matrix
mex_dir_path <- "../Input/scATAC/ME49_cell_ranger/filtered_peak_bc_matrix/"
mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '/')
feature_path <- paste(mex_dir_path, "peaks.bed", sep = '/')
barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '/')


features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)

in.dir <- "../Input/scATAC/ME49_cell_ranger/"
out.dir <- "../Input/scATAC/ME49_cell_ranger/MACS2_out/"

Tg.singlecell.info <- read.csv(paste(in.dir,"singlecell.csv", sep = ""))

mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

Tg.count.CR <- as.data.frame(as.matrix(mtx)) # Cell ranger / default
colnames(Tg.count.CR) <- gsub("\\.", "-", colnames(Tg.count.CR))
#write.csv(Pvivax.count.CR, paste(in.dir,"MACS2_out/CellrangerPvivaxPeakCounts.csv", sep = ""))


Tg.peaks.MCS <- CallPeaks(
  object = paste(in.dir, "fragments.tsv.gz", sep = ""),
  macs2.path = "/Users/kouroshz/miniconda3/envs/macs2/bin/macs2", 
  extsize = 100, 
  additional.args =" --bdg"
  
  #shift = -1,
)


# CoveragePlot(
#   object = pvivax.peaks,
#   ranges = peaks,
#   ranges.title = "MACS2"
# )

peaks.MCS <- as.data.frame(Tg.peaks.MCS)
peaks.MCS$seq.region <- 
  paste(peaks.MCS$seqnames, paste(peaks.MCS$start, peaks.MCS$end, sep = "-"), sep = "-")
peaks.MCS <- peaks.MCS %>% dplyr::select(c(seq.region, everything()))

peaks.MCS <- peaks.MCS %>% dplyr::filter(!grepl('KE', seqnames))
dim(peaks.MCS)

## peaks by cell matrix using MACS
frag.path <- paste(in.dir,'fragments.tsv.gz', sep = "")
fragments <- CreateFragmentObject(frag.path)

peak.counts <- FeatureMatrix(
  fragments = fragments,
  features = granges(Tg.peaks.MCS)
)

Tg.count.MCS <- as.data.frame(as.matrix(peak.counts))
cell.ind <- which(colnames(Tg.count.MCS) %in% colnames(Tg.count.CR))
Tg.count.MCS <- Tg.count.MCS[,cell.ind]

counts <- Tg.count.MCS %>% data.frame() %>% rownames_to_column(var = "seq.region")

sum(peaks.MCS$seq.region %in% counts$seq.region) #6110

peak_cell <- left_join(peaks.MCS, counts, by = "seq.region")

#write.csv(peak_cell, paste(out.dir,"Peaks_Cell_matrix.csv", sep = ""))


#############################################################

Tg.seqInfo <- getChromInfoFromNCBI(assembly = "TGA4", 
                                       assembled.molecules.only = TRUE, 
                                       as.Seqinfo = TRUE)

peak.cell.count <- peak_cell[,13:ncol(peak_cell)]
rownames(peak.cell.count) <- paste(peak_cell$seqnames, paste(peak_cell$start, peak_cell$end, sep = "-"),sep= ":")
colnames(peak.cell.count) <- gsub("\\.", "-", colnames(peak.cell.count))

tmp <- peak.cell.count

metadata <- read.csv(
  file = "../Input/scATAC/ME49_cell_ranger/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = tmp,
  genome = Tg.seqInfo,
  sep = c(":", "-"),
  fragments = '../Input/scATAC/ME49_cell_ranger/fragments.tsv.gz'
)

Tg <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
Tg[['peaks']]
granges(Tg)



Tg.annot <- read.table("../Input/ME49/ToxoDB-52_TgondiiME49_filter.gtf", sep = "\t")
gffRangedData <- import.gff("../Input/ME49/ToxoDB-52_TgondiiME49_filter.gtf")
annotation <- as(gffRangedData, "GRanges")
Annotation(Tg) <- annotation


# compute nucleosome signal score per cell
Tg <- NucleosomeSignal(object = Tg)
# compute TSS enrichment score per cell
Tg <- TSSEnrichment(object = Tg, tss.positions = Tg@assays$peaks@annotation,fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
Tg$pct_reads_in_peaks <- Tg$peak_region_fragments / Tg$passed_filters * 100
Tg$blacklist_ratio <- Tg$blacklist_region_fragments / Tg$peak_region_fragments

Tg$nucleosome_group <- ifelse(Tg$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = Tg, group.by = 'nucleosome_group')

VlnPlot(
  object = Tg,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


# pvivax <- SetIdent(pvivax, value = pvivax@active.ident)
# 
# pvivax <- subset(
#   x = pvivax,
#   subset = peak_region_fragments > 300 &
#     peak_region_fragments < 200 &
#     pct_reads_in_peaks > 15 &
#     blacklist_ratio < 0.05 &
#     nucleosome_signal < 4 
# )


Tg <- RunTFIDF(Tg)
Tg <- FindTopFeatures(Tg, min.cutoff = 'q0')
Tg <- RunSVD(Tg)

DepthCor(Tg)

Tg <- RunUMAP(object = Tg, reduction = 'lsi', dims = 2:30)
Tg <- FindNeighbors(object = Tg, reduction = 'lsi', dims = 2:30)
Tg <- FindClusters(object = Tg, verbose = FALSE, algorithm = 3)
DimPlot(object = Tg, label = TRUE) + NoLegend()


#gene.activities <- GeneActivity(pvivax)
gene.activities <- FeatureMatrix(
  fragments = fragments,
  features = granges(Tg.peaks.MCS)
)
peak.counts <- FeatureMatrix(
  fragments = fragments,
  features = granges(Tg.peaks.MCS)
)

Tg.count.MCS <- as.data.frame(as.matrix(peak.counts))
cell.ind <- which(colnames(Tg.count.MCS) %in% colnames(Tg.count.CR))
Tg.count.MCS <- Tg.count.MCS[,cell.ind]

Tg[['RNA']] <- CreateAssayObject(counts = Tg.count.MCS)

pvivax <- NormalizeData(
  object = pvivax,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pvivax$nCount_RNA)
)

DefaultAssay(pvivax) <- 'RNA'

FeaturePlot(
  object = pvivax,
  features = c("PvP01-14-v1-1953428-1953838","PvP01-12-v1-1940595-1941226"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)



#### Generate scHeatmap from counts
count.mat <- Tg.count.MCS
coor.dat <- data.frame(chr = unlist(lapply(strsplit(rownames(count.mat), split = '-'), `[[`, 1)),
                       strt = as.numeric(unlist(lapply(strsplit(rownames(count.mat), split = '-'), `[[`, 2))),
                       stp  = as.numeric(unlist(lapply(strsplit(rownames(count.mat), split = '-'), `[[`, 2))))


roi.chr <- "TGME49_chrIa"
roi.strt <- 928178
roi.stp <- 956595

coor.dat.slice <- coor.dat %>% dplyr::filter(chr == roi.chr & strt >= roi.strt & stp <= roi.stp)

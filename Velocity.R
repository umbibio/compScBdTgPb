library(velocyto.R)
library(Seurat)

ldat <- read.loom.matrices("../Input/compScBdTgPb/loom/27-30-33-hpi_with_10M_RH-YFP2_ME49_count.loom")
str(ldat)

S.O.intra <- readRDS('../scExpressionProfiler/S_O_toxo_MJ_labels.rds')
cell.colors <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/cell.colors.rds"))
emb <- 


hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning
# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.5)
# look at the resulting gene set
length(intersect(rownames(emat),rownames(nmat)))


library(edgeR)
library(limma)
library(openxlsx)
library(gplots)
library(dplyr)
library(tidyverse)
library(splines)

#### Test the files
#### This folder contains the output of feature count. Change as appropriate
fc.files.dir  <- "../Input/RNAseqCounts/"
### My feature coutn file names end in counts.txt, hence the grep. Change if not needed
fc.files      <- file.path(fc.files.dir, list.files(fc.files.dir)[grep("counts.txt$", list.files(fc.files.dir))])

## Creating the count table. The first column is gene ID
f <- read.table(fc.files[1], header = T, sep = '\t', stringsAsFactors = F)
f <- f[,c(1,7)]
colnames(f) <- c("GeneName", "Count1")

for(i in 2:length(fc.files)){
  tmp <- read.table(fc.files[i], header = T, sep = '\t', stringsAsFactors = F)
  tmp <- tmp[,c(1,7)]
  colnames(tmp) <- c("GeneName", paste("Count", i, sep = ''))
  f = merge(f, tmp, by.x = 1, by.y = 1)
  
}

## run5 are DNA data
RunInfo <- data.frame(Sample = gsub('_S.*', '', gsub('.counts.txt', '', basename(fc.files))), 
                      Count = paste('Count', 1:length(fc.files), sep = ''),
                      stringsAsFactors = F)
## This was manually constructed from RNAseq Info excell files
RNAseqInfo <- read.xlsx('../Input/RNAseqCounts/sampleInfo.xlsx')
RNAseqInfo <- left_join(RNAseqInfo, RunInfo, by = c('Sample.ID' = 'Sample'))

RNAseqInfo.B2 <- RNAseqInfo %>% dplyr::filter((Clone == 'B2' | Clone == 'RH' | Clone == 'B') & cond != 'fresh')
RNAseqInfo.B2$treatment <- paste(RNAseqInfo.B2$cond, RNAseqInfo.B2$passage, sep = '.')

RNAseqInfo.B2 <- RNAseqInfo.B2 %>% mutate(Names = paste(passage, cond, 'rep', Bio_rep, sep = '.'))

## Detecting Bias
y <- f %>% dplyr::select(c('GeneName', RNAseqInfo.B2$Count)) ## Just B, RH and B2 samples
## Filtering low expressions
CPM  <- cpm(y[,2:ncol(y)])
keep <-  rowSums(CPM > 2) >= 3
x <- y[keep,2:ncol(y)]
rownames(x) <- y$GeneName[keep]

treatment <- RNAseqInfo.B2$treatment[match(colnames(x), RNAseqInfo.B2$Count)]
treatment <- factor(treatment, levels = unique(treatment))

y <- DGEList(counts=x, group=treatment)
y <- calcNormFactors(y)
y$samples
#plotMDS(y)

#design <- model.matrix(~0+treatment, data=y$samples)
#colnames(design) <- levels(y$samples$group)
#y <- estimateDisp(y,design)
#plotBCV(y)



### Trying to remove batch effect from CPM values
logCPM <- cpm(y, log=TRUE, prior.count=3, normalized.lib.sizes=TRUE)
colnames(logCPM) <- RNAseqInfo.B2$Names[match(colnames(logCPM), RNAseqInfo.B2$Count)]
## MDS and heatmap clearly show a batch effect that is random
#heatmap(logCPM,cexCol = 0.5)
plotMDS(logCPM, cex = 0.5) 
## Hierarchical clustering to identify batches
clusters <- hclust(dist(t(logCPM)))
plot(clusters)
clusterCut <- cutree(clusters, 3) ## pick 3 clusters
CC <- data.frame(Samples = names(clusterCut), Batch = clusterCut, stringsAsFactors = F)
RNAseqInfo.B2 <- left_join(RNAseqInfo.B2, CC, by = c('Names' = 'Samples'))

## Remove the bad samples: These are the ones in Batch 2 and 3
bad.samples <- c("P210.extra.rep.1", "P35.extra.rep.1", "RH.intra.rep.2", "P148.extra.rep.2", 
                 "RH.extra.rep.2","P85.extra.rep.3", "P7.intra.rep.2")
RNAseqInfo.B2.bad <- RNAseqInfo.B2 %>% dplyr::filter(!(Names %in% bad.samples))


## Take good samples from Extra only
RNAseqInfo.B2.bad.extra <- RNAseqInfo.B2.bad %>% dplyr::filter(cond == 'extra' & !(passage %in% c('RH', 'P7')))

RNAseqInfo.B2.bad.extra$Time <- as.numeric(gsub("P", "", RNAseqInfo.B2.bad.extra$passage))

RNAseqInfo.B2.bad.extra <- RNAseqInfo.B2.bad.extra %>% arrange(Time)
extra.counts <- f %>% dplyr::select(RNAseqInfo.B2.bad.extra$Count)
rownames(extra.counts) <- f$GeneName

Hours <- RNAseqInfo.B2.bad.extra$Time
Time <- paste0(RNAseqInfo.B2.bad.extra$Time, "P")


y <- DGEList(counts=extra.counts, group=Time)
y$genes <- data.frame(Symbol=f$GeneName, stringsAsFactors=FALSE)

## Filter lowly expressed genes
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
plotMDS(y, labels=Hours)

## 3-knows natural splies
#X <- poly(Hours, degree=3)
#design <- model.matrix(~ X)

X <- ns(Hours, df=3)
design <- model.matrix(~X)
design

y <- estimateDisp(y, design)
sqrt(y$common.dispersion)
plotBCV(y)  

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

fit <- glmQLFTest(fit, coef=2:4)
tab <- as.data.frame(topTags(fit, n=nrow(y$counts)))
summary(decideTests(fit))

logCPM.obs <- cpm(y, log=TRUE, prior.count=fit$prior.count)
logCPM.fit <- cpm(fit, log=TRUE)


# par(mfrow=c(2,2))
# for(i in 1:4) {
#   GeneID <- row.names(tab)[i]
#   Symbol <- tab$Symbol[i]
#   logCPM.obs.i <- logCPM.obs[GeneID,]
#   logCPM.fit.i <- logCPM.fit[GeneID,]
#   plot(Hours, logCPM.obs.i, ylab="log-CPM", main=Symbol, pch=16)
#   lines(Hours, logCPM.fit.i, col="red", lwd=2)
# }

ind.sig <- tab$FDR[match(rownames(logCPM.obs), tab$Symbol)] < 0.01
print(sum(ind.sig)) ## number of DEGs

tc.logCPM <- logCPM.obs %>% as.data.frame() %>%
  mutate(GeneID = rownames(logCPM.obs)) 

tc.logCPM$trend.pval <-tab$PValue[match(tc.logCPM$GeneID, tab$Symbol)]
tc.logCPM$trend.fdr <- tab$FDR[match(tc.logCPM$GeneID, tab$Symbol)]

tc.logCPM <- tc.logCPM %>%
  pivot_longer(-c(GeneID,trend.pval, trend.fdr), names_to = "Sample", values_to = "expr")

tc.logCPM <- right_join(tc.logCPM, RNAseqInfo.B2.bad.extra, by = c('Sample' = "Count"))
tc.logCPM$rep <- tc.logCPM$Bio_rep

saveRDS(tc.logCPM, '../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_logCPM.RData')


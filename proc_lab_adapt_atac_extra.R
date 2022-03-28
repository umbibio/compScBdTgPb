library(edgeR)
library(limma)
library(openxlsx)
library(gplots)
library(dplyr)
library(tidyverse)
library(splines)

count.tab <- read.xlsx("../Input/compScBdTgPb/BulkATACToxo/Input/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_raw_counts.xlsx")
count.tab <- count.tab %>% na.omit()

count.ave <- count.tab[, 3:ncol(count.tab)]
count.ave$GeneID <- count.tab$gene_name
count.ave <- count.ave %>% group_by(GeneID) %>% summarise_all("mean")
count.ave <- count.ave %>% group_by(GeneID) %>% summarise_all("ceiling")


## Detecting Bias
y <- count.ave
## Filtering low expressions
CPM  <- cpm(y[,2:ncol(y)])
keep <-  rowSums(CPM > 2) >= 3
x <- y[keep,2:ncol(y)] %>% data.frame()
rownames(x) <- y$GeneID[keep]

treatment <- colnames(y[2:ncol(y)])
treatment <- factor(treatment, levels = unique(treatment))

y <- DGEList(counts=x, group=treatment)
y <- calcNormFactors(y)
y$samples
plotMDS(y)

extra.counts <- count.ave %>% dplyr::select(contains('P')) %>% data.frame()
rownames(extra.counts) <- count.ave$GeneID

Hours <- as.numeric(gsub("extra.P", "", colnames(extra.counts)))
Time <- paste0(Hours, "P")

ind <- sort(Hours, index.return = T)$ix

extra.counts <- extra.counts[, ind]
Hours <- Hours[ind]
Time <- Time[ind]
plotMDS(extra.counts)

y <- DGEList(counts=extra.counts, group=Time)
y$genes <- data.frame(Symbol=count.ave$GeneID, stringsAsFactors=FALSE)
plotMDS(y)

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

fit <- glmQLFit(y, design, robust=F)
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

ind.sig <- tab$FDR[match(rownames(logCPM.obs), tab$Symbol)] < 0.3
print(sum(ind.sig)) ## number of DEGs

tc.logCPM <- logCPM.obs %>% as.data.frame() %>%
  mutate(GeneID = rownames(logCPM.obs)) 

tc.logCPM$trend.pval <-tab$PValue[match(tc.logCPM$GeneID, tab$Symbol)]
tc.logCPM$trend.fdr <- tab$FDR[match(tc.logCPM$GeneID, tab$Symbol)]

tc.logCPM <- tc.logCPM %>%
  pivot_longer(-c(GeneID,trend.pval, trend.fdr), names_to = "Sample", values_to = "expr")


tc.logCPM$passage <- gsub("extra.", "", tc.logCPM$Sample)
tc.logCPM$Time <- as.numeric(gsub("extra.P", "", tc.logCPM$Sample))

saveRDS(tc.logCPM, '../Input/compScBdTgPb/LabAdaptationRNA/extra_tc_atac_logCPM.RData')


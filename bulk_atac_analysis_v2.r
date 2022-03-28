
library(openxlsx)
library(tidyverse)
library(edgeR)
library(limma)
library(splines)
library(parallel)
library(TCseq)
require(gridExtra)
library(grid)
library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(doParallel)

source("../Code/util_funcs.R")
source("../Code/util_funcsV2.R")

## Normaliization

narrow.peaks.dir <-  "~/work/sharedData/macs2_stringent/"
nr.peaks <- list.files(path = narrow.peaks.dir, pattern = ".narrowPeak")
X <- strsplit(nr.peaks, "\\.")
Names <- gsub("_peaks", "", sapply(X, "[", 1))
cond <- rep("extra", 8) 
passage <- gsub("_S[0-9]","", c(sapply(strsplit(Names, "-"), "[" , 5)[-c(7,8)], "RH","RH_gDNA"))
treatment <- paste(cond, passage, sep = ".")
Sample <- gsub(".*_", "", Names)

expmnt <- data.frame(Names, cond, passage, treatment, Sample)
expmnt
# ind <- sort(as.numeric(gsub("P", "", gsub("RH", "300", gsub("RH_gDNA", "400",expmnt$passage)))), index.return = T)$ix
# expmnt[order(expmnt$passage[ind]),]

count <- read.xlsx("../Input/BulkATACToxoPlasma/macs2_Union/peak_gene_assigned_raw_counts.xlsx")

# cpm normalization

CPM <- cbind(count[,1:2],cpm(count[,3:ncol(count)]))
CPM.avg <- CPM[,-2] %>% group_by(gene_id) %>% summarise_all("mean")
# filter low expressed genes
keep <- rowSums(CPM.avg > 2) >= 3


count.avg <- count[,-2] %>% group_by(gene_id) %>% summarise_all("mean")
count.avg <-cbind(gene_id = count.avg$gene_id, round(count.avg[,-1], digits = 0))
#count.avg <- count.avg %>% dplyr::select(-contains("gDNA"))

x <- count.avg[keep, 2:ncol(count.avg)]
rownames(x) <- count.avg$gene_id[keep]
treatment <- expmnt$treatment
treatment <- factor(treatment, levels = unique(treatment))

# Creating DGElist object

y <- DGEList(counts = x, group = treatment) 
y <- calcNormFactors(y)
y$samples

# logarithmic scale

logCPM <- cpm(y, log=TRUE, prior.count=3, normalized.lib.sizes=TRUE) 
plotMDS(logCPM, cex = 0.8)

logCPM.df <- logCPM %>% data.frame() %>% rownames_to_column(var = "GeneID")
logCPM.df <- logCPM.df %>% mutate(across(c(2:8),.fns = ~./extra.RH_gDNA))
logCPM.df <- logCPM.df %>% dplyr::select(!contains("gDNA"))
write.xlsx(logCPM.df, "../Input/BulkATACToxoPlasma/macs2_Union/Toxo_tab_ATAC_seq_logCPM.xlsx")


## Correlation and clustering
bulk.RNA <-read.xlsx('../Input/compScBdTgPb/LabAdaptationRNA/toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs_phenotype_correlations_smeV2.xlsx')
bulk.RNA.mean <- bulk.RNA  %>% dplyr::select(contains(c("GeneName",'.extra_mean')) )
bulk.RNA.mean <- bulk.RNA.mean %>% dplyr::select(!contains(c("P7")))

ind1 <- sort(as.numeric(gsub('.extra_mean', '', 
                             gsub('\\P', '', 
                                  gsub('GeneName', '1',
                                       gsub("RH",  "300", 
                                            colnames(bulk.RNA.mean)))))),index.return = T)$ix
bulk.RNA.mean <- bulk.RNA.mean[,ind1]

bulk.RNA.lng <- bulk.RNA.mean %>% pivot_longer(-GeneName, names_to = "passage", values_to = "logexpr")
bulk.RNA.lng$passage <- as.numeric(gsub("\\P", "",
                                        gsub(".extra_mean","", 
                                             gsub("RH", "300", bulk.RNA.lng$passage))))


bulk.RNA.lng <- bulk.RNA.lng %>% filter(passage != "300")  # filter out passage 300(RH)
#saveRDS(bulk.RNA.lng , "../Input/BulkATACToxoPlasma/RData/bulk_rna_seq_gene_expr_passage.RData")


bulk.ATAC <- read.xlsx("../Input/BulkATACToxoPlasma/macs2_Union/Toxo_tab_ATAC_seq_logCPM.xlsx")
bulk.ATAC <- bulk.ATAC %>% dplyr::select(!contains("P37"))

ind2 <-  sort(as.numeric(gsub("extra.", "", 
                              gsub("\\P", "", 
                                   gsub("RH", "300", 
                                        gsub("GeneID", "1", colnames(bulk.ATAC)))))), index.return = T)$ix

bulk.ATAC  <- bulk.ATAC[,ind2]
bulk.ATAC.lng <- bulk.ATAC  %>%  pivot_longer(-GeneID, names_to = "passage", values_to = "logexpr")
bulk.ATAC.lng$passage <- as.numeric(gsub("extra.RH", "300", gsub("extra.P", "", bulk.ATAC.lng$passage)))

bulk.ATAC.lng <- bulk.ATAC.lng %>% filter(passage != "300")# filter out passage 300(RH)
#saveRDS(bulk.ATAC.lng , "../Input/BulkATACToxoPlasma/RData/bulk_atac_seq_gene_expr_passage.RData")




bulk.rna <- bulk.RNA.lng
colnames(bulk.rna)[1] <- "GeneID"
bulk.atac <- bulk.ATAC.lng


comm.genes <- intersect(unique(bulk.rna$GeneID), unique(bulk.atac$GeneID))

# RNA
# no z-score scaling
bulk.rna.spline.fits <- lapply(1:length(comm.genes), function(i){
  tmp <- bulk.rna %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = passage, y = logexpr)
  bulk.rna.sp <- smooth.spline(tmp$x, tmp$y)
  bulk.rna.sp <- predict(bulk.rna.sp, seq(10, 200, by = 10)) 
  mu <- data.frame(x = bulk.rna.sp$x, y = bulk.rna.sp$y) ## For regression do not scale
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
})

bulk.rna.spline.fits <- bind_rows(bulk.rna.spline.fits)

# ATAC
# no z-score scaling
bulk.atac.spline.fits <- lapply(1:length(comm.genes), function(i){
  tmp <- bulk.atac %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = passage, y = logexpr)
  bulk.atac.sp <- smooth.spline(tmp$x, tmp$y)
  bulk.atac.sp <- predict(bulk.atac.sp, seq(10, 200, by = 10)) 
  mu <- data.frame(x = bulk.atac.sp$x, y = bulk.atac.sp$y)
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
})

bulk.atac.spline.fits <- bind_rows(bulk.atac.spline.fits)

bulk.rna.atac.spline.fits <- bind_cols(bulk.rna.spline.fits, bulk.atac.spline.fits[,2:3])
colnames(bulk.rna.atac.spline.fits) <- c('GeneID', 'x1', 'y1', 'x2', 'y2')

bulk.rna.atac.corr <- bulk.rna.atac.spline.fits %>% group_by(GeneID) %>% summarise(m.ex = mean(y1), m.ax = mean(y2), co = cor(y1, y2))

par(mfrow = c(1,1))
for(i in 1:10){
  tmp <- bulk.rna.atac.spline.fits %>% dplyr::filter(GeneID == unique(bulk.rna.atac.spline.fits$GeneID)[i])
  y.lim = c(min(c(tmp$y1, tmp$y1)), max(c(tmp$y2, tmp$y2)))
  plot(tmp$x1, tmp$y1, type = 'l', col = 'red', main = unique(bulk.rna.atac.spline.fits$GeneID)[i], ylim = y.lim)
  points(tmp$x2, tmp$y2, type = 'l', col = 'blue')
  
  # plot(tmp$x1, tmp$y1, type = 'l', col = 'red', main = paste('EXPR:', unique(sc.rna.atac.spline.fits$GeneID)[i]))
  # plot(tmp$x2, tmp$y2, type = 'l', col = 'blue', main = paste('ATAC:', unique(sc.rna.atac.spline.fits$GeneID)[i]))
  # plot(tmp$x2, tmp$y1 - tmp$y2, type = 'l', col = 'green', main = paste('Delay:', unique(sc.rna.atac.spline.fits$GeneID)[i]))
  Sys.sleep(0.9)
}


## fit spline and scale for clustering
# 
# scale.rna.spline.fits <- lapply(1:length(comm.genes), function(i){
#   tmp <- bulk.rna %>% dplyr::filter(GeneID == comm.genes[i]) %>%
#     transmute(GeneID = GeneID, x = passage, y = logexpr)
#   bulk.rna.sp <- smooth.spline(tmp$x, tmp$y)
#   bulk.rna.sp <- predict(bulk.rna.sp, seq(10, 200, by = 10))
#   mu <- data.frame(x = bulk.rna.sp$x, y = bulk.rna.sp$y) ## For regression do not scale
#   mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
#   return(mu)
# })
# 
# scale.rna.spline.fits <- bind_rows(scale.rna.spline.fits)
# 
# bulk.rna.dtw.wide <- scale.rna.spline.fits %>% 
#   pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
#   as.data.frame()
# 
# 
# scale.atac.spline.fits <- lapply(1:length(comm.genes), function(i){
#   tmp <- bulk.atac %>% dplyr::filter(GeneID == comm.genes[i]) %>%
#     transmute(GeneID = GeneID, x = passage, y = logexpr)
#   bulk.atac.sp <- smooth.spline(tmp$x, tmp$y)
#   bulk.atac.sp <- predict(bulk.atac.sp, seq(10, 200, by = 10)) 
#   mu <- data.frame(x = bulk.atac.sp$x, y = bulk.atac.sp$y) 
#   mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
#   return(mu)
# })
# 
# scale.atac.spline.fits <- bind_rows(scale.atac.spline.fits)
# 
# bulk.atac.dtw.wide <- scale.atac.spline.fits %>% 
#   pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
#   as.data.frame()

## Generate the clusters
num.clust <- 8L
bulk.rna.hc_dtws <- dtwClustCurves(bulk.rna.dtw.wide[2:ncol(bulk.rna.dtw.wide)], nclust = num.clust)
saveRDS(bulk.rna.hc_dtws, "../Input/BulkATACToxoPlasma/RData/bulk_rna_dtwclusters_8.RData")
#bulk.rna.hc_dtws <- readRDS("~/work/sharedData/bulk.rna.hc_dtws_8cluster.RData")
p <- plot(bulk.rna.hc_dtws, type = 'sc')
plot(bulk.rna.hc_dtws, type = "series", clus = 3L)
plot(bulk.rna.hc_dtws, type = "centroids", clus = 2L)


ggsave("../Input/BulkATACToxoPlasma/fig/dtw_rna_8cluster.png", 
       plot = p, 
       height = 14, width = 12, dpi = 300 )

bulk.rna.hc_dtws.df <- data.frame(GeneID = unique(colnames(bulk.rna.dtw.wide)[-1]), 
                                  order = as.numeric(bulk.rna.hc_dtws$order),
                                  cluster = cutree(bulk.rna.hc_dtws,k = num.clust))

bulk.rna.expr.hc_dtws <- left_join(bulk.rna.spline.fits, bulk.rna.hc_dtws.df, by = "GeneID")


## ggplot version of rna clusters
df.plot.rna  <- bulk.rna.expr.hc_dtws 
cluster.dtw <- paste("cluster",unique(df.plot.rna$cluster), sep = " " )

p <- lapply(unique(df.plot.rna$cluster), function(i){
  tmp <- df.plot.rna %>% dplyr::filter(cluster == i)
  ggplot(tmp, aes(x = x, y = y, color = GeneID)) + 
    geom_line(size = 0.5) + 
    theme_bw() +
    #ggtitle("rna") +
    ggtitle(cluster.dtw[i])+
    theme(strip.background = element_rect(fill = "white", colour = "black"), 
          strip.text = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.position = "none",
          axis.text = element_text(size = 16)) +
    ylab("expr") + xlab("time")
  
  
})


p1 <- do.call(grid.arrange, c(p, nrow=3))


# genes in each cluster (RNA data)

rna.clust.stat <- bulk.rna.hc_dtws.df %>% group_by(cluster) %>% summarise(total = n())


## use cluster info from rna data to identify corresponding genes in atac
bulk.rna.hc_dtws.df
bulk.atac.acc.hc_dtws <- left_join(bulk.atac.spline.fits, bulk.rna.hc_dtws.df , by = "GeneID")

atac.clust.stat <- bulk.atac.acc.hc_dtws %>% group_by(cluster) %>% distinct(GeneID) %>% summarise(total = n())

cluster.dtw <- paste("cluster",unique(bulk.atac.acc.hc_dtws$cluster), sep = " " )

df.plot.atac  <- bulk.atac.acc.hc_dtws

p <- lapply(unique(df.plot.atac$cluster), function(i){
  tmp <- df.plot.atac %>% dplyr::filter(cluster == i)
  ggplot(tmp, aes(x = x, y = y, color = GeneID)) + 
    geom_line(size = 0.5) + 
    theme_bw() +
    ggtitle(cluster.dtw[i]) +
    theme(strip.background = element_rect(fill = "white", colour = "black"), 
          strip.text = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.position = "none",
          axis.text = element_text(size = 16)) +
    ylab("expr") + xlab("time")
  
  
})

p3 <- do.call(grid.arrange, c(p, nrow=3))

p <- ggsave(filename="../Input/BulkATACToxoPlasma/fig/atac_corresponding_genes_in_dtw_rna_clusters.png",
            plot=p3,
            width = 14, height = 12,
            units = "in",
            dpi = 300)


### Regression model: From David:
##################################
# Functon to fit regression model 
# y(t) = a + b * x(t + c) + e(t) 
# by ordinary least squares
# with x and y periodic functions
##################################

rna <- bulk.rna.dtw.wide[,-1]
rownames(rna) <- bulk.rna.dtw.wide[,1]
ys <- rna

atac <- bulk.atac.dtw.wide[,-1]
rownames(atac) <- bulk.atac.dtw.wide[,1]
xs <- atac

align.xy <- function(x, y)
{
  n <- length(x)
  xbar <- mean(x)
  ybar <- mean(y)
  
  ## Center x and y
  x <- x - xbar
  y <- y - ybar
  
  ## Squared norms
  sqnrmx <- sum(x^2)
  sqnrmy <- sum(y^2)
  if (sqnrmx == 0 || sqnrmy == 0) 
    return(list(a = ybar, b = 0, c = 0, xshifted = x))
  
  ## Create shifted versions of x
  shift <- 1:n
  idx <- sapply(1:n, seq.int, length.out = n)
  xx <- c(x,x)
  xshift <- matrix(xx[idx], nrow = n)
  
  ## Estimate shift
  cp <- crossprod(y, xshift) 
  idx <- which.max(cp)
  
  ## Estimate intercept and scale
  b <- as.numeric(cp[idx]) / sqnrmx
  a <- ybar - b * xbar
  
  ## Goodness of fit
  Rsquared <- b^2 * sqnrmx / sqnrmy 
  
  list(a = a, b = b, c = idx - 1, xshifted = xshift[,idx] + xbar, 
       Rsquared = Rsquared)
}

### Test the function
# x <- sin(seq(0, pi, len=101))
# y <- 5 + 2 * x
# test <- estimate.lin.transfo(x, y)
# print(test)


## xs and ys are accessibility, expression matricies (fitted splies, with same points)
n <- ncol(xs) # number of genes
T <- nrow(xs) # number of time points
# i <- sample.int(n,1)
#matplot(tgrid, cbind(xs[,i], ys[,i]), type = "l")
out <- lapply(1:n, function(i) align.xy(xs[,i], ys[,i]))
xshift <- sapply(out, "[[", "xshifted") 
a <- sapply(out, "[[", "a")
b <- sapply(out, "[[", "b")
c <- sapply(out, "[[", "c")
# Reformulate shift in [-T/2,T/2] instead of [0,T] 
c[c > (T/2)] <- c[c > (T/2)] - T
Rsquared <- sapply(out, "[[", "Rsquared")

hist(c*6/T, main = "Distribution of shifts (Period = 6)", xlab = "Time shift")
summary(Rsquared)
png(file="../Input/BulkATACToxoPlasma/fig/atac_rna_RSquared_dist.png",
    width=600, height=350)
hist(Rsquared)
dev.off()

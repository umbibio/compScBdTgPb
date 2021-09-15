library(openxlsx)
library(tidyverse)
library(edgeR)
library(limma)
library(splines)
library(parallel)
library(sme)

## Read in microarray data from Michael White. 
count.file <-  "../Input/compScBdTgPb/syncTimeCourse/Cyclic.Genes.table.xlsx"
expmnt.file <- "../Input/compScBdTgPb/syncTimeCourse/Samples.xlsx"

count <- read.xlsx(count.file)
expmnt <- read.xlsx(expmnt.file)
GeneID <- count$GeneName
x <- count[,colnames(count) %in% expmnt$Sample[5:28]]
rownames(x) <- GeneID

## Generate time-course data
Samples <-  expmnt[5:28, ]  %>% transmute(Sample = Sample, time_point = as.numeric(gsub('h', '', Time))) %>% distinct()

## Microarray starts at S phsase. 6h corresponds to 0
Samples$time_point <- (Samples$time_point + 6) %% 12

tc_expr <- x %>% as.data.frame() %>%
  mutate(GeneID = rownames(x)) %>%
  pivot_longer(-GeneID, names_to = "Sample", values_to = "expr")

tc_expr <- right_join(tc_expr, Samples, by = 'Sample')

tc_expr.rep <- tc_expr %>% group_by(GeneID, time_point) %>% summarise(rep = 1:n())
tc_expr$rep <- tc_expr.rep$rep




saveRDS(tc_expr, '../Input/compScBdTgPb/RData/tg_tc_sync_expr.RData')


num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

sync.tc.df <- tc_expr %>% 
  transmute(y = as.numeric(expr), 
            tme = as.numeric(time_point), 
            ind = rep, variable = GeneID)


saveRDS(sync.tc.df, '../Input/compScBdTgPb/RData/tg_sync.tc.df.RData')


sync.tc.fits <- mclapply(unique(sync.tc.df$variable), 
                         function(v) 
                           sme(sync.tc.df[sync.tc.df$variable==v,c("y","tme","ind")],
                               lambda.mu = 10, lambda.v = 10), mc.cores = num.cores)

saveRDS(object = sync.tc.fits ,file = "../Input/compScBdTgPb/RData/tg_sme_fits_sync_tc_20min.RData")






## Plot a few curves to check the alignments

vs = unique(sync.tc.df$variable)[1:16]




pdf(file = "../Output/compScBdTgPb/figs/tg_sme_fits_sync.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

par(mfrow = c(4,4))
for(v in vs){
  ind <- which(unique(sync.tc.df$variable) == v)
  plot.sme(sync.tc.fits[[ind]], v, conf = F)
}

dev.off()

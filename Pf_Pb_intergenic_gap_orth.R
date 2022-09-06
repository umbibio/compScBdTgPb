library(tidyverse)
library(Biostrings)
library(bedtoolsr)
library(ggplot2)
library(openxlsx)
library(orthologr)
library(seqinr)
library(readxl)


## Read the annotations
pf.gtf <- read.table('../Input/compScBdTgPb/genes/PlasmoDB-58_Pfalciparum3D7.gtf', header = F, sep = '\t', quote = "")
pb.gtf <- read.table('../Input/compScBdTgPb/genes/PlasmoDB-58_PbergheiANKA.gtf', header = F, sep = '\t', quote = "")

all.gtf <- list(Pf = pf.gtf, Pb = pb.gtf)

all.gtf <- lapply(all.gtf, function(gtf){
  tmp <- gtf %>% dplyr::filter(V3 == 'transcript')
  tmp$GeneID = unlist(lapply(tmp$V9, function(x){
    gsub('\"', '', strsplit(x, split = ' ')[[1]][4])
  }))
  tmp <- tmp %>% transmute(chr = V1, GeneID = GeneID, strt = V4, stp = V5, strnd = V7)
  tmp <- lapply(unique(tmp$strnd), function(ss){
    tmp.ss <- tmp %>% dplyr::filter(strnd == ss) 
    tmp.ss <- tmp.ss %>% arrange(chr, strt, stp)
    tmp.ss <- lapply(unique(tmp.ss$chr), function(x){
      xx <- tmp.ss %>% dplyr::filter(chr == x)
      if(ss == '+'){
        if(nrow(xx) != 1){
          xx$upstream <- xx$strt - c(0, xx$stp[1:(nrow(xx) - 1)])
        }else{
          xx$upstream <- xx$strt - 0
        }
        
      }else if(ss == '-'){
        if(nrow(xx) != 1){
          xx$upstream <- -(xx$stp - c(xx$strt[2:nrow(xx)], (max(xx$stp) + 1000)))
        }else{
          xx$upstream <- 1000
        }
        
      }
      return(xx)
    })
    
    tmp.ss <- bind_rows(tmp.ss)
    return(tmp.ss)
  })
  
  tmp <- bind_rows(tmp)
  return(tmp)

})


pb.prod.desc <- read.csv('../Input/compScBdTgPb/genes/PBer_Prod_Desc.csv')
pb.prod.desc <- pb.prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)
all.gtf$Pb <- left_join(all.gtf$Pb, pb.prod.desc, by = 'GeneID')

pf.prod.desc <- read.csv('../Input/compScBdTgPb/genes/Pf_Prod_Desc.csv')
pf.prod.desc <- pf.prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)
all.gtf$Pf <- left_join(all.gtf$Pf, pf.prod.desc, by = 'GeneID')


write.xlsx(all.gtf$Pf, "../Input/compScBdTgPb/genes/Pf_upstream_intergenic.xlsx")
write.xlsx(all.gtf$Pb, "../Input/compScBdTgPb/genes/Pb_upstream_intergenic.xlsx")

## This requires installing Blast from NCBI

rec_Pf.vs.Pb <- blast_rec(query_file   = "../Input/compScBdTgPb/genes/PlasmoDB-58_Pfalciparum3D7_AnnotatedProteins.fasta",
                              subject_file = "../Input/compScBdTgPb/genes/PlasmoDB-58_PbergheiANKA_AnnotatedProteins.fasta",
                              delete_corrupt_cds = T, seq_type = "protein",
                              format = "fasta", blast_algorithm = "blastp",
                              eval = 0.0001, comp_cores = 16)

Pf.Pb.orth <- rec_Pf.vs.Pb %>% dplyr::select(query_id, subject_id)

Pf.Pb.orth$query_id <- gsub('\\..*', '', Pf.Pb.orth$query_id)
Pf.Pb.orth$subject_id <- gsub('\\..*', '', Pf.Pb.orth$subject_id)

colnames(Pf.Pb.orth) <- c('Pf', 'Pb')
write.xlsx(Pf.Pb.orth, "../Input/compScBdTgPb/genes/Pf_Pb_orth.xlsx")


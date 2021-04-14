library(data.table)
library(magrittr)
library(reshape2)

dt.miRNA.tissueZ.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_miRNA_log2CPM.rds')
dt.mRNA.tissueZ.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_mRNA_log2FPKM.rds')
dt.meta <- readRDS('./RDS/metadata.rds')


lst.samples.miRNA <- colnames(dt.miRNA.tissueZ.exprmat)
lst.samples.mRNA <- colnames(dt.mRNA.tissueZ.exprmat)


exprmat.mRNA.filt <- dt.mRNA.tissueZ.exprmat[,lst.samples.miRNA]
dt.meta.norm <- dt.meta[ID.sample%in%lst.samples.miRNA]



# miRNA normalization: per age Zscore -----------------------------------

dt.miRNA.ageZ.exprmat <- split(dt.meta.norm, by='Age')%>%lapply(., function(x){
  
  lst.sample.forNorm <- x$ID.sample
  exprmat.miRNA.forNorm <- dt.miRNA.tissueZ.exprmat[,lst.sample.forNorm]
  exprmat.miRNA.Norm <- apply(exprmat.miRNA.forNorm, 1, function(x){(x-mean(x))/sd(x)})
  exprmat.miRNA.Norm <- data.frame(exprmat.miRNA.Norm)
  exprmat.miRNA.Norm$ID.sample <- rownames(exprmat.miRNA.Norm)
  return(exprmat.miRNA.Norm)
  
})%>%rbindlist()

exprmat.miRNA.Zscore <- dt.miRNA.ageZ.exprmat[, 1:993]%>%as.matrix()%>%t()
colnames(exprmat.miRNA.Zscore) <- dt.miRNA.ageZ.exprmat$ID.sample



# mRNA normalization ------------------------------------------------------

dt.mRNA.ageZ.exprmat <- split(dt.meta.norm, by='Age')%>%lapply(., function(x){
  
  lst.sample.forNorm <- x$ID.sample
  exprmat.mRNA.forNorm <- dt.mRNA.tissueZ.exprmat[,lst.sample.forNorm]
  exprmat.mRNA.Norm <- apply(exprmat.mRNA.forNorm, 1, function(x){(x-mean(x))/sd(x)})
  exprmat.mRNA.Norm <- data.frame(exprmat.mRNA.Norm)
  exprmat.mRNA.Norm$ID.sample <- rownames(exprmat.mRNA.Norm)
  return(exprmat.mRNA.Norm)
  
})%>%rbindlist()

exprmat.mRNA.Zscore <- dt.mRNA.ageZ.exprmat[, 1:17249]%>%as.matrix()%>%t()
colnames(exprmat.mRNA.Zscore) <- dt.mRNA.ageZ.exprmat$ID.sample


saveRDS(exprmat.miRNA.Zscore, './RDS/exprMat/exprmat_tissue_age_Z_miRNA_log2CPM.rds')
saveRDS(exprmat.mRNA.Zscore, './RDS/exprMat/exprmat_tissue_age_Z_mRNA_log2FPKM.rds')

library(data.table)
library(magrittr)

exprmat.miRNA <- readRDS('./RDS/exprMat_miRNA_log2CPM.rds')
exprmat.mRNA <- readRDS('./RDS/exprMat_mRNA_log2FPKM.rds')
dt.meta <- readRDS('./RDS/metadata.rds')


lst.samples.miRNA <- colnames(exprmat.miRNA)
lst.samples.mRNA <- colnames(exprmat.mRNA)


exprmat.mRNA.filt <- exprmat.mRNA[,lst.samples.miRNA]
dt.meta.norm <- dt.meta[ID.sample%in%lst.samples.miRNA]




# miRNA normalization: per organ Zscore -----------------------------------

dt.miRNA.Zscore <- split(dt.meta.norm, by='Organ')%>%lapply(., function(x){
  
  lst.sample.forNorm <- x$ID.sample
  exprmat.miRNA.forNorm <- exprmat.miRNA[,lst.sample.forNorm]
  exprmat.miRNA.Norm <- apply(exprmat.miRNA.forNorm, 1, function(x){(x-mean(x))/sd(x)})
  exprmat.miRNA.Norm <- data.frame(exprmat.miRNA.Norm)
  exprmat.miRNA.Norm$ID.sample <- rownames(exprmat.miRNA.Norm)
  return(exprmat.miRNA.Norm)
 
})%>%rbindlist()

exprmat.miRNA.Zscore <- dt.miRNA.Zscore[, 1:993]%>%as.matrix()%>%t()
colnames(exprmat.miRNA.Zscore) <- dt.miRNA.Zscore$ID.sample



# mRNA normalization ------------------------------------------------------

dt.mRNA.Zscore <- split(dt.meta.norm, by='Organ')%>%lapply(., function(x){
  
  lst.sample.forNorm <- x$ID.sample
  exprmat.mRNA.forNorm <- exprmat.mRNA.filt[,lst.sample.forNorm]
  exprmat.mRNA.Norm <- apply(exprmat.mRNA.forNorm, 1, function(x){(x-mean(x))/sd(x)})%>%as.data.frame()
  exprmat.mRNA.Norm$ID.sample <- rownames(exprmat.mRNA.Norm)
  return(exprmat.mRNA.Norm)
  
})%>%rbindlist()

exprmat.mRNA.Zscore <- dt.mRNA.Zscore[, 1:17249]%>%as.matrix()%>%t()
colnames(exprmat.mRNA.Zscore) <- dt.mRNA.Zscore$ID.sample



saveRDS(exprmat.miRNA.Zscore, './RDS/exprmat_Z_miRNA_log2CPM.rds')
saveRDS(exprmat.mRNA.Zscore, './RDS/exprmat_Z_mRNA_log2FPKM.rds')

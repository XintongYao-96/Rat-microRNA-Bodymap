library(data.table)
library(magrittr)

#calculate on the non Zscore miRNA expression 
dt.miRNA.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
dt.miRNA.cor <- cor(t(dt.miRNA.exprmat))
dt.miRNA.cor.melt <- melt(dt.miRNA.cor)%>%setnames(., c('ID.Gene.A', 'ID.Gene.B', 'Corr'))%>%as.data.table()

dt.cor <- dt.miRNA.cor.melt[ID.Gene.A!=ID.Gene.B]
dt.cor[, pair := paste0(dt.cor$ID.Gene.A, '_', dt.cor$ID.Gene.B)]

saveRDS(dt.cor, './RDS/correlation_result/Corr_miRNA_to_miRNA.rds')


# calculate on the Zscore(per organ) miRNA expression 

dt.miRNA.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_miRNA_log2CPM.rds')
lst.miRNA.id <- rownames(dt.miRNA.exprmat)
lst.miRNA.ID <- gsub('\\.', '-', lst.miRNA.id)
rownames(dt.miRNA.exprmat) <- lst.miRNA.ID


dt.miRNA.cor <- cor(t(dt.miRNA.exprmat))
dt.miRNA.cor.melt <- melt(dt.miRNA.cor)%>%setnames(., c('ID.Gene.A', 'ID.Gene.B', 'Corr'))%>%as.data.table()

dt.cor <- dt.miRNA.cor.melt[ID.Gene.A!=ID.Gene.B]
dt.cor[, pair := paste0(dt.cor$ID.Gene.A, '_', dt.cor$ID.Gene.B)]

saveRDS(dt.cor, './RDS/correlation_result/Corr_miRNA_to_miRNA.Z.rds')

library(data.table)
library(magrittr)

exprmat.miRNA <- readRDS('./RDS/exprMat/exprmat_tissue_age_Z_miRNA_log2CPM.rds')
exprmat.mRNA <- readRDS('./RDS/exprMat/exprmat_tissue_age_Z_mRNA_log2FPKM.rds')

lst.sd.miRNA <- apply(exprmat.miRNA, 1, sd)
lst.miRNA <- rownames(exprmat.miRNA)[lst.sd.miRNA!=0]

lst.sd.mRNA <- apply(exprmat.mRNA, 1, sd)
lst.mRNA <- rownames(exprmat.mRNA)[lst.sd.mRNA!=0]

# Z score correlation per tissue, per age 1102-----------------------------------------------------

exprmat.all <- rbind(exprmat.miRNA, exprmat.mRNA)%>%t()
Corr.all <- cor(exprmat.all, method = 'pearson')

lst.miRNA.ID <- rownames(exprmat.miRNA)
lst.mRNA.ID <- rownames(exprmat.mRNA)

Corr.miRNA.mRNA <- Corr.all[lst.miRNA.ID, lst.mRNA.ID]

# remove mRNA included NA
lst.na.omit.mRNA <- which(colSums(is.na(Corr.miRNA.mRNA))==0)
Corr.na.omit <- Corr.miRNA.mRNA[,lst.na.omit.mRNA]


saveRDS(Corr.na.omit,'./RDS/correlation_result/correlation_tissue_age_Z.rds')

# correlation on RPM ------------------------------------------------------

exprmat.miRNA.CPM <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
exprmat.mRNA.FPKM <- readRDS('./RDS/exprMat/exprMat_mRNA_log2FPKM.rds')
exprmat.mRNA.FPKM <- exprmat.mRNA.FPKM[, colnames(exprmat.miRNA.CPM)]
lst.sd.miRNA <- apply(exprmat.miRNA.CPM, 1, sd)
lst.miRNA <- rownames(exprmat.miRNA.CPM)[lst.sd.miRNA!=0]

lst.sd.mRNA <- apply(exprmat.mRNA.FPKM, 1, sd)
lst.mRNA <- rownames(exprmat.mRNA.FPKM)[lst.sd.mRNA!=0]


exprmat.all <- rbind(exprmat.miRNA.CPM, exprmat.mRNA.FPKM)%>%t()
Corr.all <- cor(exprmat.all, method = 'pearson')


lst.miRNA.ID <- rownames(exprmat.miRNA.CPM)
lst.mRNA.ID <- rownames(exprmat.mRNA.FPKM)

Corr.miRNA.mRNA.CPM <- Corr.all[lst.miRNA.ID, lst.mRNA.ID]


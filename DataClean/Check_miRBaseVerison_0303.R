library(data.table)
library(magrittr)
dt.miRNA.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_rawFilt.rds')


dt.miRNA.ID <- fread('./tables/miRBase_miRNAList_new_old.csv')
#
identical(dt.miRNA.ID$miRNA_V22_1, dt.miRNA.ID$miRNA_V22)
identical(dt.miRNA.ID$miRNA_reference, dt.miRNA.ID$miRNA_V21)

setdiff(dt.miRNA.ID$miRNA_reference, dt.miRNA.ID$miRNA_V22)
setdiff(dt.miRNA.ID$miRNA_V22,dt.miRNA.ID$miRNA_reference)


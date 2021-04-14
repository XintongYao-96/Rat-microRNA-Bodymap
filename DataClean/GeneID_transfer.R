library(org.Rn.eg.db)
library(annotate)
library(clusterProfiler)

exprmat <- readRDS('./RDS/exprMat/exprMat_mRNA_log2FPKM.rds')
lst.genes <- rownames(exprmat)
dt.gene <- bitr(lst.genes, fromType = "ENSEMBL", #fromType是指你的数据ID类型是属于哪一类的
                toType = "SYMBOL", #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Rn.eg.db)%>%as.data.table()
saveRDS(dt.gene, './RDS/Gene_ENSEMBL_SYMBOL.rds')


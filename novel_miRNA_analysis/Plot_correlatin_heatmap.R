library(data.table)
library(magrittr)
library(pheatmap)

dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
dt.miRNA.info <- readRDS('./RDS/novel_miRNA_analysis/miRNA_info.rds')
dt.cor <- cor(dt.exprmat%>%t())
df.miRNA.info <- data.frame(dt.miRNA.info, row.names = 1)





pheatmap(dt.cor,
         show_rownames = F, show_colnames = F,
         treeheight_row = 40, treeheight_col = 40,
         annotation_col = df.miRNA.info,
         fontsize_col=12,
         filename = './charts/Novel_miRNA_correlation_heatmap.pdf'
)

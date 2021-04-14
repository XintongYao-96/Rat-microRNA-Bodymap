library(data.table)

dt.miRNA.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_miRNA_log2CPM.rds')

dt.cor <- cor(t(dt.miRNA.exprmat), method = 'pearson')
dim(dt.cor)

library(pheatmap)
library(RColorBrewer)
pheatmap(dt.cor,
         show_rownames = F, show_colnames = F,
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[8:11]))(21),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[4:8]))(99),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:4]))(21)),
         breaks = c(seq(-0.7,-0.5,0.01),
                    seq(-0.49,0.49,0.01),
                    seq(0.5,0.7,0.01)),
         border_color = NA)

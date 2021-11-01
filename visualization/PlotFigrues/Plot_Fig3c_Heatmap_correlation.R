library(data.table)
library(magrittr)
library(pheatmap)
library(RColorBrewer)

source('./scripts/PlotFigures/parameters.R')
dt.meta <- readRDS('./RDS/metadata.rds')
dt.exprmat <- readRDS('./RDS/exprMat/exprMat_logCPM_r994c318.rds')


dt.cor.sample <- cor(dt.exprmat, method = 'pearson')
annot_col <- data.frame(dt.meta[,.(ID.sample,Sex, Age,Organ)],row.names=1)
annot_color <- list(Age=colors.age,
                    Organ=colors.organ,
                    Sex = color.sex)

lst.order <- intersect(dt.meta$ID.sample, colnames(dt.cor.sample))
dt.cor.sample.order <- dt.cor.sample[lst.order, lst.order]


pheatmap(dt.cor.sample.order,
         show_rownames = FALSE, show_colnames = FALSE,
         treeheight_row = 40, treeheight_col = 40,
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         annotation_col = annot_col, annotation_colors = annot_color,
         filename = './charts/Fig3c_Heatmap_correlation.pdf', width = 8, height = 7)


source('./scripts/library.R')
source('scripts/parameters.R')

dt.cor <- readRDS('./RDS/correlation_result/correlation_tissue_Z.rds')
lst.miRNA <- gsub('.', '-', rownames(dt.cor), fixed = T)
rownames(dt.cor) <- lst.miRNA

# annotation column
lst.cluster1.mRNA <- readRDS('./RDS/correlation_result/tier1.mRNA.rds')
lst.cluster2.mRNA <- readRDS('./RDS/correlation_result/tier2.mRNA.rds')

dt.annot.column <- data.table(gene.ID = colnames(dt.cor),
                           col_type = 'none')
dt.annot.column[gene.ID%in%lst.cluster1.mRNA]$col_type <- 'cluster1'
dt.annot.column[gene.ID%in%lst.cluster2.mRNA]$col_type <- 'cluster2'
dt.annot.column$col_type <- as.factor(dt.annot.column$col_type)
df.annot.column <- data.frame(dt.annot.column, row.names = 1)


# annotation rows
lst.cluster1.miRNA <- readRDS('./RDS/correlation_result/tier1.miRNA.rds')
lst.cluster2.miRNA <- readRDS('./RDS/correlation_result/tier2.miRNA.rds')

dt.annot.row <- data.table(miRNA.ID = rownames(dt.cor),
                            row_type = 'none')
dt.annot.row[miRNA.ID%in%lst.cluster1.miRNA]$row_type <- 'cluster1'
dt.annot.row[miRNA.ID%in%lst.cluster2.miRNA]$row_type <- 'cluster2'

dt.annot.row$row_type <- as.factor(dt.annot.row$row_type)
df.annot.row <- data.frame(dt.annot.row, row.names = 1)


# color annotation
annot_col_color <- c('#07689f', '#a2d5f2', '#f0f5f9')
names(annot_col_color) <- c('cluster1', 'cluster2', 'none')

annot_row_color <- c('#e36387', '#f2aaaa', '#f0f5f9')
names(annot_row_color) <- c('cluster1', 'cluster2', 'none')

annot_color <- list(annot_col_color, annot_row_color)
names(annot_color) <- c("col_type", "row_type")


pheatmap(dt.cor,
              show_rownames = F, show_colnames = F,
              treeheight_row = 40, treeheight_col = 40,
              clustering_method = "ward.D",
              color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[8:11]))(21),
                        colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[4:8]))(99),
                        colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:4]))(21)),
              breaks = c(seq(-0.7,-0.5,0.01),
                         seq(-0.49,0.49,0.01),
                         seq(0.5,0.7,0.01)),
              annotation_row = df.annot.row,
              annotation_col = df.annot.column,
              annotation_colors = annot_color, 
              border_color = NA,
              filename = './charts/Fig6_Heatmap_withColRow.pdf', width = 18, height = 8)


# per tissue/per age Zscore -----------------------------------------------

dt.cor <- readRDS('./RDS/correlation_result/correlation_tissue_age_Z.rds')

pheatmap(dt.cor,
         show_rownames = F, show_colnames = F,
         treeheight_row = 40, treeheight_col = 40,
         clustering_method = "ward.D",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[8:11]))(21),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[4:8]))(99),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:4]))(21)),
         breaks = c(seq(-0.7,-0.5,0.01),
                    seq(-0.49,0.49,0.01),
                    seq(0.5,0.7,0.01)),
         border_color = NA,
         filename = './charts/Fig6_Heatmap_perTissueAge_Z.pdf')


# non Z score -------------------------------------------------------------

dt.cor <- readRDS('./RDS/correlation_result/correlation_nonZscore.rds')
pheatmap(t(dt.cor),
         show_rownames = F, show_colnames = F,
         treeheight_row = 40, treeheight_col = 40,
         clustering_method = "ward.D",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[8:11]))(21),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[4:8]))(99),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:4]))(21)),
         breaks = c(seq(-0.7,-0.5,0.01),
                    seq(-0.49,0.49,0.01),
                    seq(0.5,0.7,0.01)),
         border_color = NA,
         filename = './charts/Fig6_Heatmap_nonZscore.pdf')


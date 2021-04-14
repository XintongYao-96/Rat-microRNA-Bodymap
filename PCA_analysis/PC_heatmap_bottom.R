rm(list = ls())
source("./scripts/PlotFigures/parameters.R")
library(pheatmap)

dt.meta <- readRDS('./RDS/metadata.rds')
dt.met.pcs <- readRDS('./RDS/meta_PCS.rds')
dt.pca.loading <- readRDS('./RDS/PC_loadings.rds')
HCA <- readRDS('./RDS/HCA_318samples.rds')

lst.order.row <- HCA$tree_row$labels[HCA$tree_row$order]
lst.order.col <- HCA$tree_col$labels[HCA$tree_col$order]%>%intersect(.,dt.met.pcs$ID.sample)

annot_col <- data.frame(dt.meta[,.(ID.sample,Organ,Age)],row.names=1)
annot_color <- list(Organ=colors.organ,
                    Age=colors.age)


# heatmap of PC loading ---------------------------------------------------

dt.pc.mat <- dt.pca.loading[lst.order.row,]
pheatmap(dt.pc.mat[,1:5],
         show_rownames = FALSE, show_colnames = FALSE,
         scale = 'row',
         treeheight_row = 40, treeheight_col = 60,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         cluster_cols = FALSE, cluster_rows = FALSE,
         width = 2, height = 12,
         filename = './charts/Heatmap_loading.pdf' )



# heatmp of PC scores -----------------------------------------------------

dt.pc.score.mat <- dt.met.pcs[, .(ID.sample, PC1, PC2, PC3, PC4, PC5)]
Mat.pc.score.mat <- data.frame(dt.pc.score.mat, row.names = 1)%>%as.matrix(.)
Mat.pc.score.mat <- Mat.pc.score.mat[lst.order.col,]%>%t(.)

pheatmap(Mat.pc.score.mat,
         show_rownames = FALSE, show_colnames = FALSE,
         scale = 'column',
         treeheight_row = 40, treeheight_col = 60,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         cluster_cols = FALSE, cluster_rows = FALSE,
         width =12, height = 2,
         filename = './charts/Heatmap_score.pdf' )


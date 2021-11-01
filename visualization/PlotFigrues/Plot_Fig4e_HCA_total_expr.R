rm(list = ls())

library(pheatmap)
source("./scripts/PlotFigures/parameters.R")
dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.forHCA <- readRDS('./RDS/exprMat/exprMat_logCPM_r994c318.rds')



# Heatmap -----------------------------------------------------------------

IDs.miRNA <- rownames(exprMat.forHCA)
dt.miRNA.annot <- data.table(IDs.miRNA=IDs.miRNA,
                             miRNA.Type = 'Known')


dt.miRNA.annot[grep('chr', IDs.miRNA),]$miRNA.Type <- 'Novel'
annot_row <- data.frame(dt.miRNA.annot, row.names = 1)
colors.miRNA <- c('#EDA4B5', '#701473')
names(colors.miRNA) <- c('Novel', 'Known')

annot_col <- data.frame(dt.meta[,.(ID.sample,Sex,Age,Organ)],row.names=1)
annot_color <- list(Sex = color.sex,
                    Age=colors.age,
                    Organ=colors.organ,
                    miRNA.Type = colors.miRNA)

summary(as.numeric(apply(exprMat.forHCA,1,function(x)(x-mean(x))/sd(x))))

p <- pheatmap(exprMat.forHCA,
         show_rownames = FALSE, show_colnames = FALSE,
         treeheight_row = 40, treeheight_col = 60,
         clustering_method = "ward.D",scale="row",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[9:11]))(29),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[3:9]))(59),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:3]))(29)),
                  breaks = c(seq(-17,-3,0.5),
                             seq(-2.9,2.9,0.1),
                             seq(3,17,0.5)),
         annotation_row = annot_row,
                  annotation_col = annot_col, annotation_colors = annot_color,
                  annotation_names_col=FALSE,
         filename = './charts/Fig4e_HCA.pdf', width = 10 ,height = 7)


# PC loadings -------------------------------------------------------------


dt.meta <- readRDS('./RDS/metadata.rds')
dt.met.pcs <- readRDS('./RDS/meta_PCS.rds')
dt.pca.loading <- readRDS('./RDS/PC_loadings.rds')
HCA <- p

lst.order.row <- HCA$tree_row$labels[HCA$tree_row$order]
lst.order.col <- HCA$tree_col$labels[HCA$tree_col$order]%>%intersect(.,dt.met.pcs$ID.sample)

annot_col <- data.frame(dt.meta[,.(ID.sample,Organ,Age)],row.names=1)
annot_color <- list(Organ=colors.organ,
                    Age=colors.age)


# heatmap of PC loading ---------------------------------------------------

dt.pc.mat <- dt.pca.loading[lst.order.row,]
pheatmap(dt.pc.mat[,1:3],
         show_rownames = FALSE, show_colnames = FALSE,
         scale = 'row',
         treeheight_row = 40, treeheight_col = 60,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         cluster_cols = FALSE, cluster_rows = FALSE,
         width = 2, height = 12,
         filename = './charts/Fig4e_Right_PCloading.pdf' )



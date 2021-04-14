source('./scripts/parameters.R')
source('./scripts/library.R')

library(pheatmap)
dt.novel.miRNA <- readRDS('./RDS/novel_miRNA_expressed_in_organs.rds')
dt.meta <- readRDS('./RDS/metadata.rds')


exprMat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
exprMat.forHCA <- exprMat[dt.novel.miRNA$miRNA.ID, ]

annot_row <- data.frame(dt.novel.miRNA, row.names = 1)
annot_col <- data.frame(dt.meta[,.(ID.sample,Organ,Age)],row.names=1)
annot_color <- list(Organ=colors.organ,
                    Age=colors.age)

#summary(as.numeric(apply(exprMat.forHCA,1,function(x)(x-mean(x))/sd(x))))

pheatmap(exprMat.forHCA,
             show_rownames = FALSE, show_colnames = FALSE,
             treeheight_row = 40, treeheight_col = 60,
             clustering_method = "ward.D",scale="row",
             color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[9:11]))(29),
                       colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[3:9]))(59),
                       colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:3]))(29)),
             breaks = c(seq(-17,-3,0.5),
                        seq(-2.9,2.9,0.1),
                        seq(3,17,0.5)),
             annotation_col = annot_col, 
             annotation_row = annot_row,
             annotation_colors = annot_color,
             annotation_names_col=FALSE,
             filename='./charts/Fig.7.LandscapeHCA_novel.pdf',width=10,height=7)
rm(exprMat.forHCA)

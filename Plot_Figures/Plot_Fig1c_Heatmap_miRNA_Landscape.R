rm(list = ls())


source("./scripts/parameters.R")
library(pheatmap)

dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.forHCA <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')

annot_col <- data.frame(dt.meta[,.(ID.sample,Sex,Age,Organ)],row.names=1)
annot_col$Sex <- factor(annot_col$Sex, levels = c('M', 'F'))

color.cluster.type <- c('#d9ecf2','#f56a79')
names(color.cluster.type) <- c('none', 'cluster1')

annot_color <- list(Organ=colors.organ,
                    Age=colors.age,
                    type=color.cluster.type,
                    Sex=color.sex )


summary(as.numeric(apply(exprMat.forHCA,1,function(x)(x-mean(x))/sd(x))))

# heatmap without annotated miRNA (rows)
p <-pheatmap(exprMat.forHCA,
         show_rownames = FALSE, show_colnames = FALSE,
         treeheight_row = 40, treeheight_col = 60,
         clustering_method = "ward.D",scale="row",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[9:11]))(29),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[3:9]))(59),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:3]))(29)),
                  breaks = c(seq(-17,-3,0.5),
                             seq(-2.9,2.9,0.1),
                             seq(3,17,0.5)),
                  annotation_col = annot_col, annotation_colors = annot_color,
                  annotation_names_col=FALSE)
rm(exprMat.forHCA)
#saveRDS(p, './RDS/HCA_318samples.rds')


# get cluster 1 miRNA list (YXT 1028) -------------------------------------

den.miRNA <- as.dendrogram(p$tree_row)
den.cluster1.miRNA <- den.miRNA[[2]][[1]]
den.cluster2.miRNA <- den.miRNA[[2]][[2]]
order.cluster1.miRNA <- as.vector(unlist(den.cluster1.miRNA))
order.cluster2.miRNA <- as.vector(unlist(den.cluster2.miRNA))

lst.cluster1.miRNA <- p$tree_row$labels[order.cluster1.miRNA]
lst.cluster2.miRNA <- p$tree_row$labels[order.cluster2.miRNA]




annot_row <- data.table(miRNA.ID = rownames(exprMat.forHCA))
annot_row$type <- 'none'
annot_row[miRNA.ID%in%lst.cluster1.miRNA]$type <- 'cluster1'
annot_row$type <- as.factor(annot_row$type )
annot_row <- data.frame(annot_row,row.names=1)

saveRDS(lst.cluster1.miRNA, './RDS/Fig1c.cluster1.miRNA.rds')
saveRDS(lst.cluster2.miRNA, './RDS/Fig1c.cluster2.miRNA.rds')

write.table(lst.cluster1.miRNA, './supplementary/Fig1c.cluster1.miRNA.txt', quote = F, row.names = F)
write.csv(lst.cluster1.miRNA, './supplementary/SupplymentaryTable1.Fig1c.cluster1.miRNA.csv', quote = F, row.names = F)

pheatmap(exprMat.forHCA,
         show_rownames = FALSE, show_colnames = FALSE,
         treeheight_row = 20, treeheight_col = 30,
         clustering_method = "ward.D",scale="row",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[9:11]))(29),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[3:9]))(59),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:3]))(29)),
         breaks = c(seq(-17,-3,0.5),
                    seq(-2.9,2.9,0.1),
                    seq(3,17,0.5)),
         annotation_col = annot_col, annotation_colors = annot_color,
         annotation_row = annot_row,
         annotation_names_col=FALSE, 
         filename = './charts/Fig.1E.Heatmap.pdf', width = 12, height = 8)




# per tissue normalization exprmat ----------------------------------------

exprMat.forHCA <- readRDS('./RDS/exprmat_Z_miRNA_log2CPM.rds')
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
         annotation_col = annot_col, annotation_colors = annot_color,
         annotation_names_col=FALSE
         #filename='./charts/Fig.1a.LandscapeHCA.pdf',width=8,height=6
         )

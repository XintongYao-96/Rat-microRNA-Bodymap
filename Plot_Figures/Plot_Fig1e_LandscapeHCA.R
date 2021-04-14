rm(list = ls())


source("./scripts/parameters.R")
library(pheatmap)

dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.forHCA <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')

annot_col <- data.frame(dt.meta[,.(ID.sample,Sex, Age, Organ)],row.names=1)
annot_color <- list(Organ=colors.organ,
                    Age = colors.age,
                    Sex = color.sex)

summary(as.numeric(apply(exprMat.forHCA,1,function(x)(x-mean(x))/sd(x))))


# total miRNA heatmap -----------------------------------------------------


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
             annotation_names_col=FALSE,
             filename='./charts/Fig.1a.LandscapeHCA.pdf',width=10,height=6)
rm(exprMat.forHCA)
saveRDS(p, './RDS/HCA_318samples.rds')


# miRBase annotated miRNA heatmap -------------------------------------------------

# subset anntated miRNA
miRNA.info <- readRDS('./RDS/novel_miRNA_analysis/miRNA_info.rds')
lst.annotated <- miRNA.info[type.miRNA=='annotated']$miRNA.ID

exprMat.forHCA.annot <- exprMat.forHCA[lst.annotated,]

pheatmap(exprMat.forHCA.annot,
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


dt.exprmat.count <- readRDS('./RDS/exprMat/exprMat_miRNA_rawFilt.rds')

dt.exprset.count <- melt(dt.exprmat.count)%>%setnames(., c('miRNA.id', 'sample.id', 'count'))%>%as.data.table()
dt.exprset.count.filt <- dt.exprset.count[count>3]
lst.miRNA.filt <- dt.exprset.count.filt$miRNA.id%>%unique()
lst.miRNA.filt.annot <- intersect(lst.miRNA.filt, lst.annotated)

exprMat.forHCA <- exprMat.forHCA.annot[lst.miRNA.filt.annot,]

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
         annotation_names_col=FALSE)



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

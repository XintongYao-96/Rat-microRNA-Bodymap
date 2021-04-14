lst.develop.depend.miRNA <- readRDS('./RDS/list_development_dependent_miRNA.rds')

dt.meta <- readRDS('./RDS/metadata.rds')
exprMat <- readRDS('./RDS/exprMat/exprmat_Z_miRNA_log2CPM.rds')

lst.rownames <- gsub( '\\.', '-',rownames(exprMat))
rownames(exprMat) <- lst.rownames

exprMat.forHCA <- exprMat[lst.develop.depend.miRNA,]

annot_col <- data.frame(dt.meta[,.(ID.sample,Organ,Age)],row.names=1)
annot_color <- list(Organ=colors.organ,
                    Age=colors.age)


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
             annotation_names_col=FALSE,
             filename='./charts/Fig.S4.HCA_development_dependent.pdf',width=8,height=6)

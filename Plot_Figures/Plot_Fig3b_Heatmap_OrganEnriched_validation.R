rm(list = ls())


source('./scripts/parameters.R')
source('./scripts/library.R')
dt.meta <- readRDS('./RDS/metadata.rds')
#exprMat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
exprMat.Valid <- readRDS('./RDS/exprMat/exprmat_ValidationSet.rds')
dt.miRNAs.organEnriched <- readRDS('./RDS/dt.miRNAs.organEnriched')
dt.tissue.annot <- fread('./tables/Fig3_ExprMatValidate_tissue_annotation.csv')
#exprMat.valid <- fread('./tables/malerat55_miRNA_filtered_unimean_r424c180.txt')

# # make sure the miRNA name of exprmat is the same with our exprm --------
#the name rule of validation exprmat is not same with our exprmat 
#  
# lst.miRNA.rawnames <- rownames(exprMat.valid)
# lst.miRNA.names <- gsub('*', '-3p', lst.miRNA.rawnames, fixed = T)
# 
# lst.miRNA.name.5p <- paste0(lst.miRNA.names[setdiff(seq(1,424), grep('p$',lst.miRNA.names))], '-5p')
# lst.miRNA.names[setdiff(seq(1,424), grep('p$',lst.miRNA.names))] <- lst.miRNA.name.5p
# 
# lst.miRNA.ours <- rownames(exprMat)
# lst.miRNA.valid <- lst.miRNA.names
# intersect(lst.miRNA.ours, lst.miRNA.valid)
# 
# rownames(exprMat.valid) <- lst.miRNA.names
# saveRDS(exprMat.valid, './RDS/exprMat/exprmat_ValidationSet.rds')


# get miRNA list and tissue list for visualization
#IDs.miRNA.intersect <- intersect(dt.miRNAs.organEnriched$ID.gene,rownames(exprMat.Valid))
IDs.miRNA.forPlot <- c(IDs.miRNA.intersect[1:19], IDs.miRNA.intersect[26:30], IDs.miRNA.intersect[21:25])
IDs.tissue.forPlot <- dt.tissue.annot[include=='yes']$sample_name 
exprMat.Valid.forHeatmap <- exprMat.Valid[IDs.miRNA.forPlot, IDs.tissue.forPlot]


annot_row <- data.frame(dt.miRNAs.organEnriched[,.(ID.gene,Organ)],row.names = 1)

#annot_col <- data.frame(dt.tissue.annot[sample_name%in%colnames(exprMat.Valid.forHeatmap),.(sample_name, System)], row.names = 1)

annot_col <- data.frame(dt.tissue.annot[,.(sample_name, System)], row.names = 1)
annot_col$System <- factor(annot_col$System, levels=names(colors.system))
labels.col <- substr(colnames(exprMat.Valid.forHeatmap), 18, nchar(colnames(exprMat.Valid.forHeatmap)))



pheatmap(exprMat.Valid.forHeatmap,
         show_rownames = TRUE, show_colnames = TRUE,
         treeheight_row = 40, treeheight_col = 40,
         cluster_rows=FALSE, cluster_cols=FALSE,
         scale="row",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(51)),
         annotation_row = annot_row,
         annotation_col = annot_col,
         annotation_colors = list(Organ=colors.organ[unique(annot_row$Organ)],
                                  System=colors.system),
         annotation_names_col=TRUE,
         fontsize_col=12,
         labels_col = labels.col,
         border_color = F,
         filename = './charts/Fig.3b.Heatmap.Organ.Enriched.system.pdf', height = 8, width = 20)



#  combine exprmat --------------------------------------------------------

exprmat.valid <- readRDS('./RDS/exprMat/exprmat_validate_combine_miRNA_star.rds')
dt.miRNAs.organEnriched$ID.gene.group <- dt.miRNAs.organEnriched$ID.gene
dt.miRNAs.organEnriched$ID.gene.group  <- gsub('-3p', '', dt.miRNAs.organEnriched$ID.gene.group , fixed = T)
dt.miRNAs.organEnriched$ID.gene.group $ID.gene.group  <- gsub('-5p', '', dt.miRNAs.organEnriched$ID.gene.group , fixed = T)



IDs.miRNA.intersect <- intersect(lst.miRNA.organ.enrich, rownames(exprmat.valid))
IDs.miRNA.forPlot <- c(IDs.miRNA.intersect[1:6], IDs.miRNA.intersect[8:25], 
                      IDs.miRNA.intersect[31:35], IDs.miRNA.intersect[26:27], IDs.miRNA.intersect[29:30],
                      IDs.miRNA.intersect[28])
IDs.tissue.forPlot <- dt.tissue.annot[include=='yes']$sample_name 
exprMat.Valid.forHeatmap <- exprmat.valid[IDs.miRNA.forPlot, IDs.tissue.forPlot]


annot_row <- data.frame(unique(dt.miRNAs.organEnriched[,.(ID.gene.group,Organ)]),row.names = 1)
annot_col <- data.frame(dt.tissue.annot[,.(sample_name, System)], row.names = 1)
annot_col$System <- factor(annot_col$System, levels=names(colors.system))
labels.col <- substr(colnames(exprMat.Valid.forHeatmap), 18, nchar(colnames(exprMat.Valid.forHeatmap)))



pheatmap(exprMat.Valid.forHeatmap,
         show_rownames = TRUE, show_colnames = TRUE,
         treeheight_row = 40, treeheight_col = 40,
         cluster_rows=FALSE, cluster_cols=FALSE,
         scale="row",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(51)),
         annotation_row = annot_row,
         annotation_col = annot_col,
         annotation_colors = list(Organ=colors.organ[unique(annot_row$Organ)],
                                  System=colors.system),
         annotation_names_col=TRUE,
         fontsize_col=12,
         labels_col = labels.col,
         border_color = F,
         filename = './charts/Fig.3b.Heatmap.Organ.Enriched.system.pdf', height = 8, width = 20)



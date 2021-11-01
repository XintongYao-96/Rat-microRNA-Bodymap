rm(list = ls())

source('./scripts/PlotFigures/parameters.R')
source('./scripts/library.R')
library(pheatmap)
dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.Valid <- readRDS('./RDS/exprMat/exprmat_ValidationSet.rds')
dt.miRNAs.organEnriched <- readRDS('./RDS/dt.miRNAs.organEnriched.rds')
dt.tissue.annot <- fread('./tables/ExprMatValidate_tissue_annotation.csv')



# get miRNA list and tissue list for visualization ------------------------

# miRNA
IDs.miRNA.intersect <- intersect(dt.miRNAs.organEnriched$ID.gene,rownames(exprMat.Valid))
dt.miRNA.annot <- dt.miRNAs.organEnriched[ID.gene%in%IDs.miRNA.intersect,.(ID.gene,Organ)]
IDs.miRNA.forPlot <- c(dt.miRNA.annot[Organ=='Brn']$ID.gene, 
                       dt.miRNA.annot[Organ=='Tst']$ID.gene,
                       dt.miRNA.annot[Organ=='Lvr']$ID.gene,
                       dt.miRNA.annot[Organ=='Spl']$ID.gene,
                       dt.miRNA.annot[Organ=='Hrt']$ID.gene)
IDs.miRNA.forPlot <- setdiff(IDs.miRNA.forPlot, 'rno-miR-702-3p')
annot_row <- data.frame(dt.miRNA.annot, row.names = 1)

# tissue
IDs.tissue.forPlot <- dt.tissue.annot[include=='yes']$sample_name 



# Heatmap -----------------------------------------------------------------

exprMat.Valid.forHeatmap <- exprMat.Valid[IDs.miRNA.forPlot, IDs.tissue.forPlot]

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
         annotation_colors = list(Organ=colors.organ[as.character(unique(annot_row$Organ))],
                                  System=colors.system),
         annotation_names_col=TRUE,
         fontsize_col=12,
         labels_col = labels.col,
         border_color = F)
         filename = './charts/Fig5b.Heatmap_OrganEnriched_Validation.pdf', height = 8, width = 20)


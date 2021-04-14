library(data.table)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
source('./scripts/parameters.R')

dt.miRNAs.organLacked <- readRDS('./RDS/dt.miRNA.organLacked.rds')
dt.miRTarGene <- readRDS('./RDS/dt_MTI_update.rds')

# dt.MTI.IPA.Igenuity <- readRDS('./version1_scripts_results/RDS_V1/MTI/dt.MTI.IPA.Igenuity.rds')
# dt.MTI.miRecord <- readRDS('./version1_scripts_results/RDS_V1/MTI/dt.MTI.miRecord.rds')
# dt.MTI.miRTarBase <- readRDS('./version1_scripts_results/RDS_V1/MTI/dt.MTI.TarBase.rds')
# 
# lst.paies.total <- paste0(dt.miRTarGene$ID.miRNA, dt.miRTarGene$ID.gene)
# lst.pairs.IPA <- paste0(dt.MTI.IPA.Igenuity$V1, dt.MTI.IPA.Igenuity$V2)
# lst.pairs.miRrecord <- paste0(dt.MTI.miRecord$ID.miRNA, dt.MTI.miRecord$ID.gene)
# lst.pairs.miRTarbase <- paste0(dt.MTI.miRTarBase$ID.miRNA, dt.MTI.miRTarBase$ID.gene)
# 
# setdiff(lst.pairs.IPA, lst.paies.total)
# intersect(lst.pairs.miRTarbase, lst.paies.total)%>%length()



# brain -------------------------------------------------------------------

lst.brainLacked.miRNA <- dt.miRNAs.organLacked[Organ=='Brn']$ID.gene
dt.brainLacked.miRTarGene <- dt.miRTarGene[ID.miRNA%in%lst.brainLacked.miRNA]
lst.brainLacked.miRTarGene <- dt.brainLacked.miRTarGene$ID.gene


Pathway.brainLacked.miRTarGene <- data.frame(enrichGO(lst.brainLacked.miRTarGene, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))



# heatmap of target gene
exprmat.gene <- readRDS('./RDS/exprMat/exprMat_mRNA_log2FPKM.rds')

lst.gene.forHeatmap <- intersect(rownames(exprmat.gene), dt.brainLacked.miRTarGene$ID.gene)
exprmat.brainLacked.target <- exprmat.gene[lst.gene.forHeatmap,]


# heatmap
annot_col <- data.frame(dt.meta[,.(ID.sample,Organ)],row.names=1)
annot_color <- list(Organ=colors.organ)

pheatmap(exprmat.brainLacked.target,
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
         annotation_names_col=FALSE,
         filename='./charts/Fig.2b.Heatmap.BrainLacked.gene.pdf',width=8,height=6
)



#  liver ------------------------------------------------------------------

lst.liverLacked.miRNA <- dt.miRNAs.organLacked[Organ=='Lvr']$ID.gene
dt.liverLacked.miRTarGene <- dt.miRTarGene[ID.miRNA%in%lst.liverLacked.miRNA]
lst.liverLacked.miRTarGene <- dt.liverLacked.miRTarGene$ID.gene


Pathway.liverLacked.miRTarGene <- data.frame(enrichGO(lst.brainLacked.miRTarGene, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))



# heatmap of target gene

lst.gene.forHeatmap <- intersect(rownames(exprmat.gene), lst.liverLacked.miRTarGene)
exprmat.liverLacked.target <- exprmat.gene[lst.gene.forHeatmap,]


# heatmap
annot_col <- data.frame(dt.meta[,.(ID.sample,Organ)],row.names=1)
annot_color <- list(Organ=colors.organ)

pheatmap(exprmat.liverLacked.target,
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
         annotation_names_col=FALSE,
         filename='./charts/Fig.2b.Heatmap.LiverLacked.gene.pdf',width=8,height=6
)



# muscle ------------------------------------------------------------------


lst.mscLacked.miRNA <- dt.miRNAs.organLacked[Organ=='Msc']$ID.gene
dt.mscLacked.miRTarGene <- dt.miRTarGene[ID.miRNA%in%lst.mscLacked.miRNA]
lst.mscLacked.miRTarGene <- dt.mscLacked.miRTarGene$ID.gene


Pathway.mscLacked.miRTarGene <- data.frame(enrichGO(lst.brainLacked.miRTarGene, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))



# heatmap of target gene

lst.gene.forHeatmap <- intersect(rownames(exprmat.gene), lst.mscLacked.miRTarGene)
exprmat.mscLacked.target <- exprmat.gene[lst.gene.forHeatmap,]


# heatmap
annot_col <- data.frame(dt.meta[,.(ID.sample,Organ)],row.names=1)
annot_color <- list(Organ=colors.organ)

pheatmap(exprmat.mscLacked.target,
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
         annotation_names_col=FALSE,
         filename='./charts/Fig.2b.Heatmap.mscLacked.gene.pdf',width=8,height=6)
)



# utr ------------------------------------------------------------------


lst.utrLacked.miRNA <- dt.miRNAs.organLacked[Organ=='Utr']$ID.gene
dt.utrLacked.miRTarGene <- dt.miRTarGene[ID.miRNA%in%lst.utrLacked.miRNA]
lst.utrLacked.miRTarGene <- dt.utrLacked.miRTarGene$ID.gene


Pathway.utrLacked.miRTarGene <- data.frame(enrichGO(lst.utrLacked.miRTarGene, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))



# heatmap of target gene

lst.gene.forHeatmap <- intersect(rownames(exprmat.gene), lst.utrLacked.miRTarGene)
exprmat.utrLacked.target <- exprmat.gene[lst.gene.forHeatmap,]


# heatmap
annot_col <- data.frame(dt.meta[,.(ID.sample,Organ)],row.names=1)
annot_color <- list(Organ=colors.organ)

pheatmap(exprmat.utrLacked.target,
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
         annotation_names_col=FALSE,
         filename='./charts/Fig.2b.Heatmap.utrLacked.gene.pdf',width=8,height=6)




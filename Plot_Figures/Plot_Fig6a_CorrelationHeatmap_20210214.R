library(data.table)
library(magrittr)

dt.meta <- readRDS('./RDS/metadata.rds')
dt.miRNA.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_miRNA_log2CPM.rds')
dt.gene.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_mRNA_log2FPKM.rds')
lst.cluster1.miRNA <- readRDS('./RDS/correlation_result/tier1.miRNA.rds')
lst.cluster2.miRNA <- readRDS('./RDS/correlation_result/tier2.miRNA.rds')
lst.cluster1.gene <- readRDS('./RDS/correlation_result/tier1.mRNA.rds')
lst.cluster2.gene <- readRDS('./RDS/correlation_result/tier2.mRNA.rds')

rownames(dt.miRNA.exprmat) <- gsub( '.', '-',rownames(dt.miRNA.exprmat), fixed = T)


# Calculate miRNA/gene~age correlation------------------------------
# miRNA
dt.miRNA.exprmat.melt <- melt(dt.miRNA.exprmat)%>%setnames(., c('miRNA.ID','Sample.ID', 'Z_logCPM'))
dt.miRNA.exprmat.annot <- merge(dt.miRNA.exprmat.melt,dt.meta, by.x = 'Sample.ID', by.y = 'ID.sample' )%>%as.data.table()

#gene
dt.gene.exprmat.melt <- melt(dt.gene.exprmat)%>%setnames(., c('gene.ID', 'Sample.ID', 'Z_logFPKM'))
dt.gene.exprmat.annot <- merge(dt.gene.exprmat.melt,dt.meta, by.x = 'Sample.ID', by.y = 'ID.sample')%>%as.data.table()

# miRNA correlation with age
dt.miRNA.exprmat.annot$Age <- as.numeric(dt.miRNA.exprmat.annot$Age)
dt.cor.miRNA <- dt.miRNA.exprmat.annot[,  .(cor = cor(Z_logCPM, Age, method = 'spearman')), by=.(miRNA.ID)]

dt.cor.miRNA[, row_type := 'none']
dt.cor.miRNA[miRNA.ID%in%lst.cluster1.miRNA]$row_type <- 'cluster1'
dt.cor.miRNA[miRNA.ID%in%lst.cluster2.miRNA]$row_type <- 'cluster2'


#gene correlation with age
dt.gene.exprmat.annot$Age <- as.numeric(dt.gene.exprmat.annot$Age)
dt.cor.gene <- dt.gene.exprmat.annot[,  .(cor = cor(Z_logFPKM, Age, method = 'spearman')), by=.(gene.ID)]

dt.cor.gene[, col_type := 'none']
dt.cor.gene[gene.ID%in%lst.cluster1.gene]$col_type <- 'cluster1'
dt.cor.gene[gene.ID%in%lst.cluster2.gene]$col_type <- 'cluster2'




# Plot heatmap ------------------------------------------------------------

source('./scripts/library.R')
source('scripts/parameters.R')

dt.cor <- readRDS('./RDS/correlation_result/correlation_tissue_Z.rds')
lst.miRNA <- gsub('.', '-', rownames(dt.cor), fixed = T)
rownames(dt.cor) <- lst.miRNA

# annotation column
dt.annot.column <- dt.cor.gene[gene.ID%in%colnames(dt.cor)]
dt.annot.column$col_type <- as.factor(dt.annot.column$col_type)
df.annot.column <- data.frame(dt.annot.column, row.names = 1)


# annotation rows

dt.annot.row <- dt.cor.miRNA[miRNA.ID%in%rownames(dt.cor)]
dt.annot.row$row_type <- as.factor(dt.annot.row$row_type)
df.annot.row <- data.frame(dt.annot.row, row.names = 1)


# color annotation
annot_col_color <- c('#07689f', '#a2d5f2', '#f0f5f9')
names(annot_col_color) <- c('cluster1', 'cluster2', 'none')

annot_row_color <- c('#e36387', '#f2aaaa', '#f0f5f9')
names(annot_row_color) <- c('cluster1', 'cluster2', 'none')

annot_cor_color <- brewer.pal(11, 'PiYG')

annot_color <- list(col_type=annot_col_color,
                    row_type=annot_row_color,
                    cor = annot_cor_color)

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
         filename = './charts/Fig6_Heatmap_withAgeCor.pdf', width = 15, height = 8)




# Plot_Fig7b_correlation_violin -------------------------------------------

dt.cor.forPlot <- 


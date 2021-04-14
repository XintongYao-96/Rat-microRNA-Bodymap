rm(list = ls())
library(data.table)
source('./scripts/library.R')
source('./scripts/parameters.R')
dt.DEGs.forw <- readRDS('./RDS/DEG_results/DEG_OrganPerAge_org.rds')

dt.DEGs.rev <- data.table(
  logFC = -dt.DEGs.forw$logFC,
  AveExpr = dt.DEGs.forw$AveExpr,
  t = dt.DEGs.forw$t,
  P.Value = dt.DEGs.forw$P.Value,
  adj.P.Val = dt.DEGs.forw$adj.P.Val,
  B = dt.DEGs.forw$B,
  ID.gene = dt.DEGs.forw$ID.gene,
  Age = dt.DEGs.forw$Age,
  Organ.A = dt.DEGs.forw$Organ.B,
  Organ.B = dt.DEGs.forw$Organ.A
)

dt.DEGs <- rbindlist(list(dt.DEGs.forw,dt.DEGs.rev))

rm(dt.DEGs.forw,dt.DEGs.rev)


# Identify organ-lacked miRNAs --------------------------------------------
cutOff=1
dt.miRNAs.organLacked <- dt.DEGs[adj.P.Val<0.05][logFC<(-log2(cutOff))][,.N,by=.(Organ.B,ID.gene)][N==40]
setnames(dt.miRNAs.organLacked,'Organ.B','Organ')

saveRDS(dt.miRNAs.organLacked[,.(ID.gene,Organ)],'./RDS/dt.miRNA.organLacked.rds')


# Draw heatmap ------------------------------------------------------------

# Set row orders

organs.Lacked.orderByLackN <- dt.miRNAs.organLacked[,.N,by='Organ'][order(-N,Organ)]$Organ%>%as.character()
organs.orderByLackN <- c(organs.Lacked.orderByLackN,setdiff(organs,organs.Lacked.orderByLackN))

dt.miRNAs.organLacked$Organ <- factor(dt.miRNAs.organLacked$Organ,levels=organs.Lacked.orderByLackN)
ID.miRNAs.OrganLacked.orderByOrgan <- dt.miRNAs.organLacked[order(Organ,ID.gene)]$ID.gene

# Set col orders
dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.miRNA.logCPM <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')


IDs.sample.forPlot <- colnames(exprMat.miRNA.logCPM)
dt.meta$Organ <- factor(dt.meta$Organ,levels=organs.orderByLackN)

IDs.sample.forPlot.orderByOrgan <- dt.meta[IDs.sample.forPlot][order(Organ,Age,Sex,BioRep)]$ID.sample


# annotation bar
annot_row <- data.frame(dt.miRNAs.organLacked[,.(ID.gene,Organ)],row.names = 1)
annot_col <- data.frame(dt.meta[,.(ID.sample,Organ)],row.names=1)

# plot
library(pheatmap)
library(RColorBrewer)


pheatmap(exprMat.miRNA.logCPM[ID.miRNAs.OrganLacked.orderByOrgan,IDs.sample.forPlot.orderByOrgan],
         show_rownames = TRUE, show_colnames = FALSE,
         treeheight_row = 40, treeheight_col = 40,
         cluster_rows=FALSE, cluster_cols=FALSE,
         scale="row",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[9:11]))(29),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[3:9]))(59),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:3]))(29)),
         breaks = c(seq(-18,-3,0.5),
                    seq(-2.9,2.9,0.1),
                    seq(3,18,0.5)),
         annotation_col = annot_col,
         annotation_row = annot_row,
         annotation_colors = list(Organ=colors.organ),
         annotation_names_col=FALSE,annotation_names_row=FALSE,
         fontsize_row=12,
        filename = './charts/Fig.2a.Heatmap.Organ-lacked.pdf',width=8,height=5)  

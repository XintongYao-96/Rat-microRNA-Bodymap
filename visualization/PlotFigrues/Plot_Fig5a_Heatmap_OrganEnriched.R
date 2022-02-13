rm(list = ls())

source('./scripts/PlotFigures/parameters.R')


# define organ enriched miRNAs --------------------------------------------

dt.DEGresult <- readRDS('./RDS/DEG_OrganPerAge_combine.rds')

# find organ enriched miRNAs
cutOff=4
dt.miRNAs.organEnriched.FC4 <- dt.DEGresult[adj.P.Val<0.05][logFC>(log2(cutOff))][,.N,by=.(Organ.B,ID.gene)][N==40]
setnames(dt.miRNAs.organEnriched.FC4,'Organ.B','Organ')
# saveRDS(dt.miRNAs.organEnriched.FC4 , './RDS/dt.miRNAs.organEnriched.rds')

dt.miRNAs.organEnriched <- dt.miRNAs.organEnriched.FC4
# saveRDS(dt.miRNAs.organEnriched,'./RDS/dt.miRNAs.organEnriched.rds')



# Plot --------------------------------------------------------------------


dt.meta <- readRDS('./RDS/metadata.rds')

# Set row orders
organs.Enriched.orderByN <- dt.miRNAs.organEnriched[,.N,by='Organ'][order(-N,Organ)]$Organ
organs.all.orderByN <- c(organs.Enriched.orderByN,setdiff(organs,organs.Enriched.orderByN))

dt.miRNAs.organEnriched$Organ <- factor(dt.miRNAs.organEnriched$Organ,levels=organs.Enriched.orderByN)
ID.miRNAs.organEnriched.orderByOrgan <- dt.miRNAs.organEnriched[order(Organ,ID.gene)]$ID.gene

# Set col orders
exprMat.miRNA.logCPM <- readRDS('./RDS/exprMat/exprMat_logCPM_r604c318.rds')
IDs.sample.forPlot <- colnames(exprMat.miRNA.logCPM)
dt.meta$Organ <- factor(dt.meta$Organ,levels=organs.all.orderByN)
IDs.sample.forPlot.orderByOrgan <- dt.meta[IDs.sample.forPlot][order(Organ,Age,Sex,BioRep)]$ID.sample

# annotation bar
annot_row <- data.frame(dt.miRNAs.organEnriched[,.(ID.gene,Organ)],row.names = 1)
annot_col <- data.frame(dt.meta[,.(ID.sample,Organ)],row.names=1)

# subset expression mat
exprMat.forPlot <- exprMat.miRNA.logCPM[ID.miRNAs.organEnriched.orderByOrgan ,IDs.sample.forPlot.orderByOrgan]
exprMat.forPlot.zscore <- apply(exprMat.forPlot,1,function(x)(x-mean(x))/sd(x))

bound <- ceiling(max(abs(as.numeric(exprMat.forPlot.zscore))))
core <- min(bound,3)

pheatmap(exprMat.forPlot,
         show_rownames = TRUE, show_colnames = FALSE,
         treeheight_row = 40, treeheight_col = 40,
         cluster_rows=FALSE, cluster_cols=FALSE,
         scale="row", # per miRNA score
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[9:11]))((bound-core)/0.1+1),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[3:9]))(2*core/0.1-1),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:3]))((bound-core)/0.1+1)
         ),
         breaks = c(seq(-bound,-core,0.1),
                    seq(-core+0.1,core-0.1,0.1),
                    seq(core,bound,0.1)),
         annotation_col = annot_col,
         annotation_row = annot_row,
         annotation_colors = list(Organ=colors.organ),
         annotation_names_col=FALSE,annotation_names_row=FALSE,
         fontsize_row=12,
         border_color = NA,
         filename = './charts/Fig5a_Heatmap_OrganEnriched.pdf',width=20,height=12)  

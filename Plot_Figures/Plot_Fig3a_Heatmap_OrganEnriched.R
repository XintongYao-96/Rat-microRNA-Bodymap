rm(list = ls())

source('./scripts/parameters.R')
source('./scripts/library.R')


# define organ enriched miRNAs --------------------------------------------

dt.DEGresult <- readRDS('./RDS/DEG_results/DEG_OrganPerAge_combine.rds')

# find organ enriched miRNAs
cutOff=4
dt.miRNAs.organEnriched.FC4 <- dt.DEGresult[adj.P.Val<0.05][logFC>(log2(cutOff))][,.N,by=.(Organ.B,ID.gene)][N==40]
setnames(dt.miRNAs.organEnriched.FC4,'Organ.B','Organ')
saveRDS(dt.miRNAs.organEnriched.FC4 , './RDS/dt.miRNAs.organEnriched.rds')

dt.miRNAs.organEnriched <- dt.miRNAs.organEnriched.FC4


# Plot --------------------------------------------------------------------


dt.meta <- readRDS('./RDS/metadata.rds')

# Set row orders
organs.Enriched.orderByN <- dt.miRNAs.organEnriched[,.N,by='Organ'][order(-N,Organ)]$Organ
organs.all.orderByN <- c(organs.Enriched.orderByN,setdiff(organs,organs.Enriched.orderByN))

dt.miRNAs.organEnriched$Organ <- factor(dt.miRNAs.organEnriched$Organ,levels=organs.Enriched.orderByN)
ID.miRNAs.organEnriched.orderByOrgan <- dt.miRNAs.organEnriched[order(Organ,ID.gene)]$ID.gene

# Set col orders
exprMat.miRNA.logCPM <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')


IDs.sample.forPlot <- colnames(exprMat.miRNA.logCPM)
dt.meta$Organ <- factor(dt.meta$Organ,levels=organs.all.orderByN)

IDs.sample.forPlot.orderByOrgan <- dt.meta[IDs.sample.forPlot][order(Organ,Age,Sex,BioRep)]$ID.sample

# annotation bar
annot_row <- data.frame(dt.miRNAs.organEnriched[,.(ID.gene,Organ)],row.names = 1)
annot_col <- data.frame(dt.meta[,.(ID.sample,Organ)],row.names=1)

# subset expression mat
exprMat.forPlot <- exprMat.miRNA.logCPM[ID.miRNAs.organEnriched.orderByOrgan ,IDs.sample.forPlot.orderByOrgan]

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
         filename = './charts/Fig.3a.Heatmap.OrganEnriched.pdf',width=20,height=12)  


# check the expression level of miR-124 （2021.02.25）-----------------------------------

lst.samples.brn <- dt.meta[Organ=='Brn']$ID.sample
exprmat.brn <- exprMat.forPlot[, intersect(colnames(exprMat.forPlot), lst.samples.brn)]
lst.mean.expr <- apply(exprmat.brn, 1, mean)
sort(lst.mean.expr)


dt.exprSet <- melt(dt.exprmat)
setnames(dt.exprSet, c('miRNA.ID', 'Sample.ID', 'logCPM'))
dt.meta.exprSet <- merge(dt.exprSet, dt.meta, by.x = 'Sample.ID', by.y = 'ID.sample')%>%as.data.table()
lst.miRNA.let7 <- grep('rno-let-7*', dt.meta.exprSet$miRNA.ID, value = T)%>%unique()
dt.exprSet.let7 <- dt.meta.exprSet[miRNA.ID%in%lst.miRNA.let7]


ggplot(dt.exprSet.let7, aes(x=Organ, y =logCPM))+
  geom_violin(trim = F)

ggplot(dt.exprSet.let7, aes(x=Age, y =logCPM))+
  geom_violin(trim = F)

ggplot(dt.exprSet.let7, aes(x=Sex, y =logCPM))+
  geom_violin(trim = F)

ggplot(dt.exprSet.let7, aes(x=miRNA.ID, y =logCPM))+
  geom_violin(trim = F)+
  geom_jitter(aes(color = Organ))+
  theme(axis.text.x = element_text(size=12,face='bold',color='gray40', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, face = 'bold'),
        axis.title = element_text(size=18,face='bold'),
        legend.title = element_text(size=18,face='bold'),
        legend.text = element_text(size=18,face='bold'))
  

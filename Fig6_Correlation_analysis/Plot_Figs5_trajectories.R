library(data.table)
library(magrittr)
library(ggplot2)

lst.cluster1.miRNA <- readRDS('./RDS/correlation_result/tier1.miRNA.rds')
lst.cluster2.miRNA <- readRDS('./RDS/correlation_result/tier2.miRNA.rds')
lst.cluster1.gene <- readRDS('./RDS/correlation_result/tier1.mRNA.rds')
lst.cluster2.gene <- readRDS('./RDS/correlation_result/tier2.mRNA.rds')
dt.meta <- readRDS('./RDS/metadata.rds')
dt.miRNA.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_miRNA_log2CPM.rds')
dt.gene.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_mRNA_log2FPKM.rds')

rownames(dt.miRNA.exprmat) <- gsub( '.', '-',rownames(dt.miRNA.exprmat), fixed = T)



# annot expression level with metadata ------------------------------
# miRNA
dt.miRNA.exprmat.melt <- melt(dt.miRNA.exprmat)%>%setnames(., c('miRNA.ID','Sample.ID', 'Z_logCPM'))
dt.miRNA.exprmat.annot <- merge(dt.miRNA.exprmat.melt,dt.meta, by.x = 'Sample.ID', by.y = 'ID.sample' )%>%as.data.table()


#gene
dt.gene.exprmat.melt <- melt(dt.gene.exprmat)%>%setnames(., c('gene.ID', 'Sample.ID', 'Z_logFPKM'))
dt.gene.exprmat.annot <- merge(dt.gene.exprmat.melt,dt.meta, by.x = 'Sample.ID', by.y = 'ID.sample')%>%as.data.table()


# get cluster 1 & 2 expression --------------------------------------------

dt.cluster1.miRNA.exprmat.annot <- dt.miRNA.exprmat.annot[miRNA.ID%in%lst.cluster1.miRNA]
dt.cluster2.miRNA.exprmat.annot <- dt.miRNA.exprmat.annot[miRNA.ID%in%lst.cluster2.miRNA]

dt.cluster1.miRNA.exprmat.annot$group <- paste0(dt.cluster1.miRNA.exprmat.annot$Organ, 
                                                dt.cluster1.miRNA.exprmat.annot$miRNA.ID, 
                                                dt.cluster1.miRNA.exprmat.annot$BioRep, 
                                                dt.cluster1.miRNA.exprmat.annot$Sex)

dt.cluster2.miRNA.exprmat.annot$group <- paste0(dt.cluster2.miRNA.exprmat.annot$Organ, 
                                                dt.cluster2.miRNA.exprmat.annot$miRNA.ID, 
                                                dt.cluster2.miRNA.exprmat.annot$BioRep, 
                                                dt.cluster2.miRNA.exprmat.annot$Sex)
# dt.cluster1.miRNA.exprmat.annot$development <- as.character(dt.cluster1.miRNA.exprmat.annot$Age)%>%as.numeric(.)
# dt.cluster2.miRNA.exprmat.annot$Age <- as.numeric(dt.cluster2.miRNA.exprmat.annot$Age)


# cluster 1 miRNA expression
dt.cluster1.miRNA.exprmat.annot$Age <- as.numeric(dt.cluster1.miRNA.exprmat.annot$Age)
p <- ggplot(dt.cluster1.miRNA.exprmat.annot, aes(x= Age, y=Z_logCPM, group = group))+
  geom_line(color = 'gray80')+
  labs(x='Age', y= 'Zscore of logCPM')+
  theme_bw()+
  labs(x = 'Age (Week)', y = 'Z score of miRNA expression \n across 11 organs and both sexes',
       title = 'miRNA trajectories of cluster 1')+
  theme(plot.title = element_text(size=20,face='bold',color='gray40', hjust = 0.5),
        axis.text.x = element_text(size=12,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 12, face = 'bold',    color='gray40'),
        axis.title = element_text(size=12,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
ggsave('./charts/Fig.S6.cluster1.miRNA.expression.pdf', p, width = 6, height = 5)

dt.cluster1.miRNA.exprmat.annot$Age <- as.numeric(dt.cluster1.miRNA.exprmat.annot$Age)
lst.cluster1.miRNA.age.cor <- cor(dt.cluster1.miRNA.exprmat.annot$Z_logCPM, dt.cluster1.miRNA.exprmat.annot$Age, method = 'spearman')

dt.cluster2.miRNA.exprmat.annot$Age <- as.numeric(dt.cluster2.miRNA.exprmat.annot$Age)
lst.cluster2.miRNA.age.cor <- cor(dt.cluster2.miRNA.exprmat.annot$Z_logCPM, dt.cluster2.miRNA.exprmat.annot$Age, method = 'spearman')

dt.cluster2.miRNA.exprmat.annot$Organ <- as.factor(dt.cluster2.miRNA.exprmat.annot$Organ)
dt.cluster2.miRNA.exprmat.annot$Sex <- as.factor(dt.cluster2.miRNA.exprmat.annot$Sex)
dt <- dt.cluster2.miRNA.exprmat.annot[,  .(cor = cor(Z_logCPM, Age)), by=.(Organ, Sex)]




# cluster 2 miRNA expression
ggplot(dt.cluster2.miRNA.exprmat.annot, aes(x= Age, y=Z_logCPM, group = group))+
  geom_line(color = 'gray80')+
  labs(x='Age', y= 'Zscore of logCPM')+
  theme_bw()+
  labs(x = 'Age (Week)', y = 'Z score of miRNA expression \n across 11 organs and both sexes',
       title = 'miRNA trajectories of cluster 2')+
  theme(plot.title = element_text(size=20,face='bold',color='gray40', hjust = 0.5),
        axis.text.x = element_text(size=12,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 12, face = 'bold',    color='gray40'),
        axis.title = element_text(size=12,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
ggsave('./charts/Fig.S6.cluster2.miRNA.expression.pdf', width = 6, height = 5)


#  gene -------------------------------------------------------------------

dt.gene.exprmat.annot$group <- paste0(dt.gene.exprmat.annot$Organ, 
                                      dt.gene.exprmat.annot$gene.ID, 
                                      dt.gene.exprmat.annot$BioRep, 
                                      dt.gene.exprmat.annot$Sex)

dt.cluster1.gene.exprmat.annot <- dt.gene.exprmat.annot[gene.ID%in%lst.cluster1.gene]
dt.cluster2.gene.exprmat.annot <- dt.gene.exprmat.annot[gene.ID%in%lst.cluster2.gene]

# cluster 1 gene
p <- ggplot(dt.cluster1.gene.exprmat.annot, aes(x= Age, y=Z_logFPKM, group = group))+
  geom_line(color = 'gray80')+
  labs(x='Age', y= 'Zscore of logCPM')+
  theme_bw()+
  labs(x = 'Age (Week)', y = 'Z score of miRNA expression \n across 11 organs and both sexes',
       title = 'Gene trajectories of cluster 1')+
  theme(plot.title = element_text(size=20,face='bold',color='gray40', hjust = 0.5),
        axis.text.x = element_text(size=12,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 12, face = 'bold',    color='gray40'),
        axis.title = element_text(size=12,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
ggsave('./charts/Fig.S6c.cluster1.gene.pdf', p, width = 6, height = 5)

# cluster 2 gene
p <- ggplot(dt.cluster2.gene.exprmat.annot, aes(x= Age, y=Z_logFPKM, group = group))+
  geom_line(color = 'gray80')+
  labs(x='Age', y= 'Zscore of logCPM')+
  theme_bw()+
  labs(x = 'Age (Week)', y = 'Z score of miRNA expression \n across 11 organs and both sexes',
       title = 'Gene trajectories of cluster 2')+
  theme(plot.title = element_text(size=20,face='bold',color='gray40', hjust = 0.5),
        axis.text.x = element_text(size=12,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 12, face = 'bold',    color='gray40'),
        axis.title = element_text(size=12,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
ggsave('./charts/Fig.S6d.cluster2.gene.pdf', p, width = 6, height = 5)


dt.cluster1.gene.exprmat.annot$Age <- as.numeric(dt.cluster1.gene.exprmat.annot$Age)
lst.cluster1.gene.age.cor <- cor(dt.cluster1.gene.exprmat.annot$Z_logFPKM, dt.cluster1.gene.exprmat.annot$Age, method = 'spearman')

dt.cluster2.gene.exprmat.annot$Age <- as.numeric(dt.cluster2.gene.exprmat.annot$Age)
lst.cluster2.gene.age.cor <- cor(dt.cluster2.gene.exprmat.annot$Z_logFPKM, dt.cluster2.gene.exprmat.annot$Age, method = 'spearman')




# gene expression in interest pathway
dt.pathway.gene <- readRDS('./RDS/Pathway_result/dt.pathway.gene.rds')
dt.pathway.cluster1.gene <- readRDS('./RDS/Pathway_result/dt.pathway.cluster1.gene.rds')
dt.pathway.cluster2.gene <- readRDS('./RDS/Pathway_result/dt.pathway.cluster2.gene.rds')


dt.pathway.cluster1.gene.annot <- dt.gene.exprmat.annot[gene.ID%in%dt.pathway.cluster1.gene$ID.gene]

# Figure s6e: pathway 1 gene 
p <- ggplot(dt.pathway.cluster1.gene.annot, aes(x= Age, y=Z_logFPKM, group = group))+
  geom_line(color = 'gray80')+
  labs(x='Age', y= 'Zscore of logCPM')+
  theme_bw()+
  labs(x = 'Age (Week)', y = 'Z score of miRNA expression \n across 11 organs and both sexes',
       title = 'Trajectories of cell division relevant genes')+
  theme(plot.title = element_text(size=18,face='bold',color='gray40', hjust = 0.5),
        axis.text.x = element_text(size=12,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 12, face = 'bold',    color='gray40'),
        axis.title = element_text(size=12,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
ggsave('./charts/FigS6e.pathway1.gene.expression.pdf', width = 6, height = 5)



# Figure s6f: pathway 2 gene 

dt.pathway.cluster2.gene.annot <- dt.gene.exprmat.annot[gene.ID%in%dt.pathway.cluster2.gene$ID.gene]

p <- ggplot(dt.pathway.cluster2.gene.annot, aes(x= Age, y=Z_logFPKM, group = group))+
  geom_line(color = 'gray80')+
  labs(x='Age', y= 'Zscore of logCPM')+
  theme_bw()+
  labs(x = 'Age (Week)', y = 'Z score of miRNA expression \n across 11 organs and both sexes',
       title = 'Trajectories of immune response relevant genes')+
  theme(plot.title = element_text(size=16,face='bold',color='gray40', hjust = 0.5),
        axis.text.x = element_text(size=12,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 12, face = 'bold',    color='gray40'),
        axis.title = element_text(size=12,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
ggsave('./charts/FigS6f.pathway2.gene.expression.pdf', p, width = 6, height = 5)


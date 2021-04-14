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


# miRNA correlation with age
dt.miRNA.exprmat.annot$Age <- as.numeric(dt.miRNA.exprmat.annot$Age)
dt.cluster1.miRNA.exprmat.annot <- dt.miRNA.exprmat.annot[miRNA.ID%in%lst.cluster1.miRNA]
dt.cluster2.miRNA.exprmat.annot <- dt.miRNA.exprmat.annot[miRNA.ID%in%lst.cluster2.miRNA]

dt.cor.clu1.miRNA <- dt.cluster1.miRNA.exprmat.annot[,  .(cor = cor(Z_logCPM, Age, method = 'spearman')), by=.(Organ, Sex)]
dt.cor.clu1.miRNA$cluster <- 'miRNA.Cluster1'

dt.cor.clu2.miRNA <- dt.cluster2.miRNA.exprmat.annot[,  .(cor = cor(Z_logCPM, Age, method = 'spearman')), by=.(Organ, Sex)]
dt.cor.clu2.miRNA$cluster <- 'miRNA.Cluster2'


# gene correlation with age 
dt.gene.exprmat.annot$Age <- as.numeric(dt.gene.exprmat.annot$Age)
dt.cluster1.gene.exprmat.annot <- dt.gene.exprmat.annot[gene.ID%in%lst.cluster1.gene]
dt.cluster2.gene.exprmat.annot <- dt.gene.exprmat.annot[gene.ID%in%lst.cluster2.gene]

dt.cor.clu1.gene <- dt.cluster1.gene.exprmat.annot[,  .(cor = cor(Z_logFPKM, Age)), by=.(Organ, Sex)]
dt.cor.clu1.gene$cluster <- 'Gene.Cluster1'

dt.cor.clu2.gene <- dt.cluster2.gene.exprmat.annot[,  .(cor = cor(Z_logFPKM, Age)), by=.(Organ, Sex)]
dt.cor.clu2.gene$cluster <- 'Gene.Cluster2'


# bind correlation and plot
dt.cor <- rbind(dt.cor.clu1.miRNA,
                dt.cor.clu2.miRNA,
                dt.cor.clu1.gene,
                dt.cor.clu2.gene)
dt.cor$cluster <- factor(dt.cor$cluster, levels = c('miRNA.Cluster1', 'miRNA.Cluster2', 'Gene.Cluster1', 'Gene.Cluster2'))

lst.color <- c( '#e36387', '#f2aaaa','#07689f', '#a2d5f2')

ggplot(dt.cor, aes(x=cluster, y=cor))+
  geom_violin(trim = F, aes(fill=cluster), color = 'white')+
  geom_boxplot(width=.1)+
  scale_fill_manual(values = lst.color)+
  theme_bw()+
  labs(x='Cluster', y= 'Spearman correlation')+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=20,face='bold', hjust = 1, angle = 45),
        axis.text.y = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size=20,face='bold'),
        legend.title = element_text(size=20,face='bold'),
        legend.text = element_text(size=20,face='bold'))
ggsave('./charts/Fig.6.violin.cor.withAge.pdf', width = 6, height = 10 )  

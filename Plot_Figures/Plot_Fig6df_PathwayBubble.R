library(data.table)
library(magrittr)
library(ggplot2)
library(magrittr)
library(RColorBrewer)

dt.cluster1.pathway <- readRDS('./RDS/Pathway_result/Pathway.tier1.gene.rds')
dt.cluster2.pathway <- readRDS('./RDS/Pathway_result/Pathway.tier2.gene.rds')


# pathway of interest (10 cluster1 pathway+18 cluster 2 pathway) ----------

dt.cluster1.top10.pathway <- data.table(dt.cluster1.pathway)[1:10,]
dt.cluster1.top10.pathway$cluster <- 'cluster1'
dt.cluster2.top10.pathway <- data.table(dt.cluster2.pathway)[1:10,]
dt.cluster2.top10.pathway$cluster <- 'cluster2'
dt.cluster2.top9.pathway <- dt.cluster2.top10.pathway[Description != 'adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains']
  
dt.cluster.pathway <- rbind(dt.cluster1.top10.pathway, dt.cluster2.top9.pathway)



# generate Fold enrichment ----------------------------------------------

numerator.GeneRatio <- strsplit(dt.cluster.pathway$GeneRatio, '/')%>%sapply(., '[',1)%>%as.numeric()
denominator.GeneRatio <- strsplit(dt.cluster.pathway$GeneRatio, '/')%>%sapply(., '[',2)%>%as.numeric()
numberator.BgRatio <- strsplit(dt.cluster.pathway$BgRatio, '/')%>%sapply(., '[',1)%>%as.numeric()
denominator.BgRatio <- strsplit(dt.cluster.pathway$BgRatio, '/')%>%sapply(., '[',2)%>%as.numeric()
dt.cluster.pathway$FC.enrich <- (numerator.GeneRatio/denominator.GeneRatio)/(numberator.BgRatio/denominator.BgRatio)
dt.cluster.pathway$logP <- -log10(dt.cluster.pathway$p.adjust)
dt.cluster.pathway$Description <- as.factor(dt.cluster.pathway$Description)
dt.cluster.pathway.forPlot <- dt.cluster.pathway[, .(Description, Count, cluster, FC.enrich, logP)]

# bubble plot -------------------------------------------------------------

color.low<- brewer.pal(11, 'RdYlBu')[2]
color.mid <- brewer.pal(11, 'RdYlBu')[6]

color.high <- brewer.pal(11, 'RdYlBu')[10]

ggplot(dt.cluster.pathway[cluster == 'cluster1'], aes(x=FC.enrich, y=Description))+
  geom_point(aes(size=Count, colour = logP))+
  scale_color_gradient(low = '#f1d4d4',
                       high = '#6a097d')+
  theme_bw()+
  scale_x_continuous(limits = c(2,5))+
  labs(x='Fold Enrichment', 
       y= 'Pathway')+
  theme(axis.text.x = element_text(size=12,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 12, face = 'bold',    color='gray40'),
        axis.title = element_text(size=12,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
ggsave('./charts/Fig.6b.PathwayBubble.cluster1gene.pdf', width = 6, height = 5)


ggplot(dt.cluster.pathway.forPlot[cluster == 'cluster2'], aes(x=FC.enrich, y=Description))+
  geom_point(aes(size=Count, colour = logP))+
  scale_color_gradient(low = '#f1d4d4',
                       high = '#6a097d')+
  theme_bw()+
  labs(x='Fold Enrichment', 
       y= 'Pathway')+
  scale_x_continuous(limits = c(2,9))+
  theme(axis.text.x = element_text(size=12,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 12, face = 'bold',    color=s'gray40'),
        axis.title = element_text(size=12,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
ggsave('./charts/Fig.6b.PathwayBubble.cluster2gene.pdf', width = 8, height = 5)

dt.cluster.pathway[Description=='antigen processing and presentation']


library(data.table)
library(magrittr)
source('./scripts/parameters.R')
dt.miRNA.pattern <- readRDS('./RDS/miRNA_development_pattern_proportion.rds')
dt.miRNA.pattern$omic <- 'miRNA'
dt.miRNA.pattern$Perc <-dt.miRNA.pattern$Perc
dt.gene.pattern <- fread('./tables/RNAseq_development pattern.csv')

dt.gene.pattern.melt <- melt(dt.gene.pattern)%>%setnames(., c('Organ', 'Pattern', 'Perc'))
dt.gene.pattern.melt$omic <- 'gene'

dt.pattern <- rbind(dt.miRNA.pattern, dt.gene.pattern.melt)
dt.pattern$omic <- factor(dt.pattern$omic, levels = c('miRNA', 'gene'))
library(ggplot2)

ggplot(dt.pattern, aes(x=omic, y=Perc))+
  geom_boxplot()+
  facet_wrap(~Pattern, nrow=3)+
  geom_jitter(aes(colour=Organ),width = 0.1)+
  scale_color_manual(values = colors.organ)+
  theme_bw()

dt.pattern.forPlot <- dt.pattern[Pattern%in%c('MMM', 'DMM','UMM','DDM','MDM')]
dt.pattern.forPlot$Pattern <- factor(dt.pattern.forPlot$Pattern, levels = c('MMM', 'DMM','UMM','DDM','MDM'))

# figure 3c
ggplot(dt.pattern.forPlot, aes(x=omic, y=Perc))+
  geom_boxplot()+
  facet_wrap(~Pattern, nrow=1)+
  geom_jitter(aes(colour=Organ),width = 0.1)+
  scale_color_manual(values = colors.organ)+
  theme_bw()+
  labs(x= 'Dataset', y= 'Percentage of genes/miRNAs 
       within each pattern')+
  theme(axis.text.x = element_text(size=16,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 16, face = 'bold',    color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=16,face='bold'),
        legend.text = element_text(size=16,face='bold',color='gray40'),
        strip.text = element_text(size=16,face='bold'),
        legend.position = 'none')
ggsave('./charts/Figure4c.miRNA.gene.pattern.percentage.pdf', width = 15, height = 4)

dt.pattern.forPlot[Pattern == 'MMM'][omic == 'miRNA']$Perc%>%mean()
dt.pattern.forPlot[Pattern == 'MMM'][omic == 'gene']$Perc%>%mean()



# supplementary figure 4
ggplot(dt.pattern, aes(x=omic, y=Perc))+
  geom_boxplot()+
  facet_wrap(~Pattern, nrow=3)+
  geom_jitter(aes(colour=Organ),width = 0.1)+
  scale_color_manual(values = colors.organ)+
  theme_bw()+
  labs(x= 'Dataset', y= 'Percentage of genes/miRNAs 
       within each pattern')+
  theme(axis.text.x = element_text(size=16,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 16, face = 'bold',    color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=16,face='bold'),
        legend.text = element_text(size=16,face='bold',color='gray40'),
        strip.text = element_text(size=16,face='bold'),
        legend.position = 'none')
ggsave('./charts/FigureS4.miRNA.gene.pattern.percentage.pdf', width = 25, height = 10)

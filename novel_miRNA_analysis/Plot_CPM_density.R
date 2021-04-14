library(magrittr)
library(data.table)
library(ggplot2)


dt.miRNA.info <- readRDS('./RDS/novel_miRNA_analysis/miRNA_info.rds')
setkey(dt.miRNA.info, 'miRNA.ID')
exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_CPM.rds')



exprSet <- melt(exprmat)%>%setnames(., c('miRNA.ID', 'Sample.ID', 'logCPM'))
dt.exprSet <- data.table(exprSet, key = 'miRNA.ID' )
dt.exprSet.info <- dt.exprSet[dt.miRNA.info]

ggplot(dt.exprSet.info, aes(x=logCPM, group = type.miRNA))+
  geom_density(aes(color = type.miRNA, fill = type.miRNA), alpha = 0.5)+
  scale_x_log10()+
  theme_bw()+
  labs(x= 'log10 (CPM)', y = 'Density')+
  theme(axis.text.x = element_text(size=20,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 20, face = 'bold',    color='gray40'),
        axis.title = element_text(size=20,face='bold'),
        legend.title = element_text(size=20,face='bold'),
        legend.text = element_text(size=20,face='bold',color='gray40'))
ggsave('./charts/Novel_log10CPM_density.pdf', width = 8, height = 6)

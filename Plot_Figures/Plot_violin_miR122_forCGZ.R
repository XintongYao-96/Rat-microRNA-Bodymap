library(data.table)
library(magrittr)
library(ggplot2)
source('./scripts/library.R')

dt.meta <- readRDS('./RDS/metadata.rds')
dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')

# expression and annotation
dt.exprset <- melt(dt.exprmat)%>%setnames(., c('miRNA.ID', 'Sample.ID', 'logCPM'))%>%as.data.table()
dt.annot.exprset <- merge(dt.meta, dt.exprset,
                          by.x = "ID.sample", by.y = 'Sample.ID',
                          all = F)

# generate order by 122 expression level
dt.annot.exprset[miRNA.ID=='rno-miR-122-5p', mean.miR122 := mean(logCPM), by = .(Organ)]
dt.annot.exprset.miR122 <- dt.annot.exprset[miRNA.ID=='rno-miR-122-5p', .(Organ, mean.miR122)]%>%unique()
level.organ <- dt.annot.exprset.miR122[order(-mean.miR122)]$Organ
dt.annot.exprset$Organ <- factor(dt.annot.exprset$Organ , level.organ)

ggplot(dt.annot.exprset[miRNA.ID=='rno-miR-122-5p'], aes(x= Organ, y = logCPM))+
  geom_violin(trim = F, width = 1.5, aes(fill = Organ,color = Organ))+
  geom_boxplot(width = .1, color = 'gray40')+
  scale_color_manual(values = colors.organ)+
  scale_fill_manual(values = colors.organ)+
  labs(title = 'Expression of miR-122 across 11 organs of rats')+
  theme_bw()+  
  theme(legend.position = 'none',
    plot.title = element_text(size = 18, face = 'bold', hjust = 0.5),
    axis.text.x = element_text(size=18,face='bold',color='gray40', hjust = 0.5),
    axis.text.y = element_text(size = 18, face = 'bold'),
    axis.title = element_text(size=18,face='bold'),
    legend.title = element_text(size=18,face='bold'),
    legend.text = element_text(size=18,face='bold'))
ggsave('./charts/Violin_miR122_expression_11organs.pdf', width = 8, height = 5)
               
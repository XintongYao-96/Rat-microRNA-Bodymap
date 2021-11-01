rm(list = ls())

library(data.table)
library(magrittr)
library(ggplot2)

dt.meta <- readRDS('./RDS/metadata.rds')
dt.exprSet <- readRDS('./RDS/exprMat/exprSet_rawcount_1033miR320Sample.rds')
dt.miRNA.info <- readRDS('./RDS/miRNA_Chromosome_annot.rds')


dt.exprSet.novel <- dt.exprSet[miRNA.ID%in%dt.miRNA.info[type.miRNA=='novel']$miRNA.name]

dt.exprSet.novel.meta <-merge(dt.exprSet.novel, dt.meta, 
                              by.x = 'Sample.ID', by.y =  "Colnames.miRNA",
                              all = T)

dt.exprSet.expressed <- dt.exprSet.novel.meta[Count>3]
length(unique(dt.exprSet.expressed$miRNA.ID))


dt.expressed.times <- dt.exprSet.expressed[, .(N.expressed.times = .N), by = .(miRNA.ID, Organ)]
# dt.expressed.times[miRNA.ID=='rno-chr7-4593']
dt.expressed.organ.times <- dt.expressed.times[N.expressed.times>3, .(N.organs = .N), by = miRNA.ID]
# saveRDS(dt.expressed.organ.times, './RDS/novel_miRNA_expressed_in_organs.rds')


lst.novel.nonExpressed <- setdiff(dt.exprSet.novel$miRNA.ID%>%unique(), unique(dt.exprSet.expressed$miRNA.ID))
dt.exprmat <- readRDS('./RDS/exprMat/exprMat_rawcount_r994c318.rds')
dt.exprmat[lst.novel.nonExpressed,]%>%apply(., 1, sum)

dt.forPlot.novel <- dt.expressed.organ.times[, .(N = .N), by = N.organs]
ggplot(dt.forPlot.novel, aes(x= N.organs, y = N))+
  geom_bar(stat = 'identity', fill = '#9a8194')+
  scale_x_continuous(breaks = seq(1, 11, 1))+
  labs(x = 'Number of organs', y = 'Number of miRNAs')+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(size=18,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 18, face = 'bold'),
        axis.title = element_text(size=18,face='bold'),
        legend.title = element_text(size=18,face='bold'),
        legend.text = element_text(size=18,face='bold'))
# ggsave('./charts/Barplot_number_of_novel_miRNA_detected_in_organs.pdf', width = 4,height = 8)  



# annotated miRNAs  -------------------------------------------------------

dt.exprSet.novel <- dt.exprSet[miRNA.ID%in%dt.miRNA.info[type.miRNA=='novel']$miRNA.name]

dt.exprSet.novel.meta <-merge(dt.exprSet.novel, dt.meta, 
                              by.x = 'Sample.ID', by.y =  "Colnames.miRNA",
                              all = T)


dt.exprSet.annot <- dt.exprSet[miRNA.ID%in%dt.miRNA.info[type.miRNA=='annotated']$miRNA.name]

dt.exprSet.annot.meta <-merge(dt.exprSet.annot, dt.meta, 
                              by.x = 'Sample.ID', by.y = "Colnames.miRNA",
                              all= T)
dt.exprSet.expressed <- dt.exprSet.annot.meta[Count>3]
length(unique(dt.exprSet.expressed$miRNA.ID))


dt.expressed.times <- dt.exprSet.expressed[, .(N.expressed.times = .N), by = .(miRNA.ID, Organ)]
# dt.expressed.times[miRNA.ID=='rno-chr7-4593']
dt.expressed.organ.times <- dt.expressed.times[N.expressed.times>3, .(N.organs = .N), by = miRNA.ID]
dt.forPlot.annot <- dt.expressed.organ.times[, .(N = .N), by = N.organs]
# 
# ggplot(dt.forPlot.annot, aes(x= N.organs, y = N))+
#   geom_bar(stat = 'identity', fill = '#9a8194')+
#   scale_x_continuous(breaks = seq(1, 11, 1))+
#   labs(x = 'Number of organs', y = 'Number of miRNAs')+
#   coord_flip()+
#   theme_bw()+
#   theme(axis.text.x = element_text(size=18,face='bold',color='gray40', hjust = 0.5),
#         axis.text.y = element_text(size = 18, face = 'bold'),
#         axis.title = element_text(size=18,face='bold'),
#         legend.title = element_text(size=18,face='bold'),
#         legend.text = element_text(size=18,face='bold'))
# 



# Plot annot and novel together -------------------------------------------

dt.forPlot.annot[, perc := N/sum(N)]
dt.forPlot.annot[, value := -perc*100]
dt.forPlot.annot[, type := 'annotated']
dt.forPlot.novel[, perc := N/sum(N)]
dt.forPlot.novel[, value := perc*100]
dt.forPlot.novel[, type := 'novel']

dt.forPlot <- rbind(dt.forPlot.annot, dt.forPlot.novel)

ggplot(dt.forPlot, aes(x= N.organs, y = value))+
  geom_bar(stat = 'identity',pos='stack',color='gray80',
           aes(fill=type))+
  scale_x_continuous(breaks = seq(1, 11, 1))+
  labs(x = 'Number of organs', y = 'Proportion of miRNAs detected (%)')+
  # scale_y_continuous(limits = c(-55, 40))+
  scale_fill_manual(values = rev(c('#a2d5f2','#07689f')))+
  # coord_flip()+
  theme_bw()+
  theme(
        axis.text.x = element_text(size=20,face='bold',color='black', hjust = 0.5),
        axis.text.y = element_text(size = 20, color='black', face = 'bold'),
        axis.title = element_text(size=20,color='black', face='bold'),
        legend.title = element_text(size=20,face='bold'),
        legend.text = element_text(size=20,face='bold'))
ggsave('./charts/Fig7d_WaterfulPlot_miRNA_organ.pdf', width = 15, height = 6)


# stats -------------------------------------------------------------------

dt.forPlot[N.organs<4][type=='novel']$perc%>%sum()
dt.forPlot[N.organs<4][type=='annotated']$perc%>%sum()

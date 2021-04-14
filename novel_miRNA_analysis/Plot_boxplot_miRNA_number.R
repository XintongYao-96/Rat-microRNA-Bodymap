######################################
# Compare miRNA NUMBER               #
# NOVEL V.S. ANNOTATED               #
#  Xintong                           #
######################################




library(data.table)
library(magrittr)
source('./scripts/parameters.R')
dt.meta <- readRDS('./RDS/metadata.rds')
dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_rawFilt.rds')
dt.miRNA.info <- readRDS('./RDS/novel_miRNA_analysis/miRNA_info.rds')

# generate novel miRNA expression
lst.novel.miRNA <- dt.miRNA.info[type.miRNA=='novel']$miRNA.ID
dt.novel.exprmat <- dt.exprmat[lst.novel.miRNA, ]

dt.exprset.novel <- melt(dt.novel.exprmat)%>%setnames(., c('miRNA.ID', 'Sample.ID', 'count'))%>%as.data.table()
dt.exprset.expressed.novel <- dt.exprset.novel[count>3]
saveRDS(dt.exprset.novel, './RDS/exprMat/exprSet_novel_count.rds')

# generate annot miRNA expression
lst.annot.miRNA <- dt.miRNA.info[type.miRNA=='annotated']$miRNA.ID
dt.annot.exprmat <- dt.exprmat[lst.annot.miRNA, ]

# count >3 are detected miRNAs
dt.exprset.annot <- melt(dt.annot.exprmat)%>%setnames(., c('miRNA.ID', 'Sample.ID', 'count'))%>%as.data.table()
dt.exprset.expressed.annot <- dt.exprset.annot[count>3]
saveRDS(dt.exprset.annot, './RDS/exprMat/exprSet_annot_count.rds')


dt.N.miRNA.novel <- dt.exprset.expressed.novel[, .(N.miRNA.novel = .N), by = 'Sample.ID']
dt.N.miRNA.annot <- dt.exprset.expressed.annot[, .(N.miRNA.annot = .N), by = 'Sample.ID']

dt.N.miRNA <- merge(dt.N.miRNA.novel, dt.N.miRNA.annot, by = 'Sample.ID')%>%as.data.table()
dt.meta.N.miRNA <- merge(dt.meta, dt.N.miRNA, by.x = 'ID.sample', by.y = 'Sample.ID')
dt.meta.N.miRNA[, Perc.Novel := N.miRNA.novel/(N.miRNA.novel+N.miRNA.annot)]



# Scatterplot 0207 correlation of annot miRNA to novel miRNA
ggplot(dt.meta.N.miRNA, aes(x = N.miRNA.annot, y= N.miRNA.novel))+
  geom_point(aes(color = Organ),size = 3)+
  scale_color_manual(values = colors.organ)+
  geom_smooth(method = 'lm')+
  labs(x = 'Number of miRNA annotated in miRBase',
       y = 'Number of novel miRNA ')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=14,face='bold',color='gray40', vjust = 0.5),
        axis.text.y = element_text(size = 14, face = 'bold',    color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=20,face='bold'),
        legend.text = element_text(size=20,face='bold',color='gray40'))
ggsave('./charts/Novel_to_annot_scatterplot.pdf', width = 5, height = 5)  



# organize the boxplot in an order of perc.novel
dt.median <- dt.meta.N.miRNA[, .(median(Perc.Novel)), by = Organ]%>%setorder(., - V1)
order.organ <- dt.median$Organ
dt.meta.N.miRNA$Organ <- factor(dt.meta.N.miRNA$Organ, levels = order.organ)

ggplot(dt.meta.N.miRNA, aes(x = Organ, y= Perc.Novel))+
  geom_boxplot(aes(fill = Organ))+
  scale_color_manual(values = colors.organ)+
  scale_fill_manual(values = colors.organ)+
  theme_bw()+
  labs(y = 'Percentage of novel miRNA',
       x = 'Organ')+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=14,face='bold',color='gray40', vjust = 0.5),
        axis.text.y = element_text(size = 14, face = 'bold',    color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=18,face='bold'),
        legend.text = element_text(size=18,face='bold',color='gray40'))
ggsave('./charts/Novel_percent_boxplot.pdf', width = 8, height = 5)  




# generate dtforPlot

dt.meta.N.miRNA$Group <- paste0(dt.meta.N.miRNA$Organ, dt.meta.N.miRNA$Age)

dt.meta.N.miRNA.mean <- dt.meta.N.miRNA[, .(mean = mean(N.miRNA.novel)), by= c('Group')]
dt.meta.N.miRNA.sd <- dt.meta.N.miRNA[, .(sd = sd(N.miRNA.novel)), by= c('Group')]



dt.meta.N.miRNA.forPlot <- merge(dt.meta.N.miRNA.mean, dt.meta.N.miRNA.sd, by = 'Group')
dt.meta.N.miRNA.forPlot$Organ <- factor(substr(dt.meta.N.miRNA.forPlot$Group, 1,3))
dt.meta.N.miRNA.forPlot$Age <- factor(substr(dt.meta.N.miRNA.forPlot$Group, 4, nchar(dt.meta.N.miRNA.forPlot$Group)),
                                      levels = c('2','6','21','104'))

levels.group <- paste0(rep(names(colors.organ), each = 4), rep(c('2', '6', '21', '104'),11))
dt.meta.N.miRNA.forPlot$Group <- factor(dt.meta.N.miRNA.forPlot$Group, levels = levels.group)


ggplot(dt.meta.N.miRNA.forPlot, aes(x= Group, y= mean, fill=Organ))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(x= Group, ymin=mean-sd, ymax=mean+sd),width=0.1, color='gray40', position = position_dodge(0.9),  # 设置误差线颜色，宽度等
                size=0.8)+
  labs(x = '11 organs across 4 ages', y = 'Number of miRNA detected')+
  facet_grid(~ Organ, scales = 'free_x', space = 'free_x')+
  scale_fill_manual(values = colors.organ)+
  theme(panel.border = element_blank(),
        strip.text.x = element_text(size=0),
        panel.spacing.x = unit(0.2, "mm")
        # strip.background = element_rect(color='black')
  )+
  theme(axis.text.x = element_text(size=14,face='bold',color='gray40', angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 14, face = 'bold',    color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=18,face='bold'),
        legend.text = element_text(size=18,face='bold',color='gray40'))

ggsave('./charts/Novel.BarPlot.N.miRNA.pdf', width = 15, height = 8)



dt.meta.N.miRNA[Organ=='Lvr']$Perc.Novel%>%mean()





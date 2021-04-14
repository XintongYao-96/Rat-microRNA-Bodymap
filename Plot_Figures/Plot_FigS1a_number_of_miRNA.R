library(data.table)
library(magrittr)
source('./scripts/parameters.R')
dt.meta <- readRDS('./RDS/metadata.rds')
dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_rawFilt.rds')
hist(dt.exprmat)

dt.exprset <- melt(dt.exprmat)%>%setnames(., c('miRNA.ID', 'Sample.ID', 'count'))%>%as.data.table()
dt.exprset.expressed <- dt.exprset[count>3]

dt.N.miRNA <- dt.exprset.expressed[, .(N.miRNA = .N), by = 'Sample.ID']
dt.meta.N.miRNA <- merge(dt.meta, dt.N.miRNA, by.x = 'ID.sample', by.y = 'Sample.ID')
dt.meta.N.miRNA$Group <- paste0(dt.meta.N.miRNA$Organ, dt.meta.N.miRNA$Age)

dt.meta.N.miRNA.mean <- dt.meta.N.miRNA[, .(mean = mean(N.miRNA)), by= c('Group')]
dt.meta.N.miRNA.sd <- dt.meta.N.miRNA[, .(sd = sd(N.miRNA)), by= c('Group')]


# generate dtforPlot
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

ggsave('./charts/FigS1a.BarPlot.N.miRNA.pdf', width = 15, height = 8)




# statistic for text ------------------------------------------------------


mean(dt.meta.N.miRNA.forPlot$mean)
mean(dt.meta.N.miRNA.forPlot[Organ=='Lvr']$mean)

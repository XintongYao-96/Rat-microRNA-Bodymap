library(data.table)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)

dt.meta <- readRDS('./RDS/metadata.rds')
dt.exprmat <- readRDS('./RDS/exprMat/exprMat_logCPM_r994c318.rds')

dt.cor.sample <- cor(dt.exprmat, method = 'pearson')
dt.melt.cor <- melt(dt.cor.sample)%>%setnames(., c('ID.SampleA', 'ID.SampleB', 'cor'))%>%as.data.table()
dt.melt.cor <- dt.melt.cor[ID.SampleA!=ID.SampleB]


# specify comparing type --------------------------------------------------

dt.melt.cor$Biological.GroupA <- substr(dt.melt.cor$ID.SampleA, 1, 9)
dt.melt.cor$Biological.GroupB <- substr(dt.melt.cor$ID.SampleB, 1, 9)

dt.melt.cor$Organ_ageA <- paste0(substr(dt.melt.cor$Biological.GroupA, 1,3),substr(dt.melt.cor$Biological.GroupA, 7,9))
dt.melt.cor$Organ_ageB <- paste0(substr(dt.melt.cor$Biological.GroupB, 1,3),substr(dt.melt.cor$Biological.GroupB, 7,9))

dt.melt.cor$OrganA <- substr(dt.melt.cor$Biological.GroupA, 1,3)
dt.melt.cor$OrganB <- substr(dt.melt.cor$Biological.GroupB, 1,3)

dt.melt.cor$Biological.Group.Type <- ifelse(dt.melt.cor$Biological.GroupA==dt.melt.cor$Biological.GroupB , 'Biological_Replicate', 
                                           ifelse(dt.melt.cor$Organ_ageA==dt.melt.cor$Organ_ageB, 'Same_organ_age',
                                                  ifelse(dt.melt.cor$OrganA==dt.melt.cor$OrganB, 'Same_organ_cross_age', 'Different_organ')))

dt.melt.cor$Biological.Group.Type <- factor(dt.melt.cor$Biological.Group.Type, levels =
                                             c("Biological_Replicate", "Same_organ_age", "Same_organ_cross_age", "Different_organ"))

dt.melt.cor$Biological.Group.Type <- ifelse(dt.melt.cor$Biological.GroupA==dt.melt.cor$Biological.GroupB , 'Biological_Replicate', 'Different_Sample')

ggplot(dt.melt.cor, aes(x=Biological.Group.Type, y= cor))+
  geom_violin(aes(fill = Biological.Group.Type, color = Biological.Group.Type))+
  geom_boxplot(color = 'gray40', width = .05)+
  labs(x= 'Group Type', 
       y = 'Pearson correlation coefficient')+
  theme_bw()+
  scale_fill_manual(values = c('#B1DBEC', '#E07594'))+
  scale_color_manual(values = c('#B1DBEC', '#E07594'))+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=15, face='bold',color='black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15, face = 'bold',  color='black'),
        axis.title = element_text(size=15,face='bold'),
        legend.title = element_text(size=16,face='bold'),
        legend.text = element_text(size=15,face='bold',color='gray40'))
ggsave('./charts/Fig3a_Boxplot_correlation.pdf', width = 3.6, height = 6)




# stat --------------------------------------------------------------------

mean.rep <- mean(dt.melt.cor[Biological.Group.Type == 'Biological_Replicate']$cor)
sd.rep <- sd(dt.melt.cor[Biological.Group.Type == 'Biological_Replicate']$cor)


mean.diff <- mean(dt.melt.cor[Biological.Group.Type == 'Different_Sample']$cor)
sd.rep <- sd(dt.melt.cor[Biological.Group.Type == 'Different_Sample']$cor)

t.test(dt.melt.cor[Biological.Group.Type == 'Biological_Replicate']$cor,
       dt.melt.cor[Biological.Group.Type == 'Different_Sample']$cor)

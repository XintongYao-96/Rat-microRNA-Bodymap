library(data.table)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)

dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')

dt.cor.sample <- cor(dt.exprmat, method = 'spearman')

dt.meta <- readRDS('./RDS/metadata.rds')

annot_col <- data.frame(dt.meta[,.(ID.sample,Organ,Age)],row.names=1)
annot_color <- list(Organ=colors.organ,
                    Age=colors.age)

pheatmap(dt.cor.sample,
         show_rownames = FALSE, show_colnames = FALSE,
         treeheight_row = 40, treeheight_col = 40,
         clustering_method = "ward.D",
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         annotation_col = annot_col, annotation_colors = annot_color,
         filename = './charts/FigS1_Heatmap_correlation.pdf')
      
         






dt.melt.cor <- melt(dt.cor.sample)%>%setnames(., c('ID.SampleA', 'ID.SampleB', 'cor'))%>%as.data.table()
dt.melt.cor <- dt.melt.cor[ID.SampleA!=ID.SampleB]


dt.melt.cor$Biological.GroupA <- substr(dt.melt.cor$ID.SampleA, 1, 9)
dt.melt.cor$Biological.GroupB <- substr(dt.melt.cor$ID.SampleB, 1, 9)

dt.melt.cor$Organ_ageA <- paste0(substr(dt.melt.cor$Biological.GroupA, 1,3),substr(dt.melt.cor$Biological.GroupA, 7,9))
dt.melt.cor$Organ_ageB <- paste0(substr(dt.melt.cor$Biological.GroupB, 1,3),substr(dt.melt.cor$Biological.GroupB, 7,9))

dt.melt.cor$OrganA <- substr(dt.melt.cor$Biological.GroupA, 1,3)
dt.melt.cor$OrganB <- substr(dt.melt.cor$Biological.GroupB, 1,3)

#dt.melt.cor$Biological.Group.Type <- ifelse(dt.melt.cor$Biological.GroupA==dt.melt.cor$Biological.GroupB , 'Biological_Replicate', 
#                                            ifelse(dt.melt.cor$Organ_ageA==dt.melt.cor$Organ_ageB, 'Same_organ_age',
#                                                   ifelse(dt.melt.cor$OrganA==dt.melt.cor$OrganB, 'Same_organ_cross_age', 'Different_organ')))

#dt.melt.cor$Biological.Group.Type <- factor(dt.melt.cor$Biological.Group.Type, levels = 
#                                              c("Biological_Replicate", "Same_organ_age", "Same_organ_cross_age", "Different_organ"))

dt.melt.cor$Biological.Group.Type <- ifelse(dt.melt.cor$Biological.GroupA==dt.melt.cor$Biological.GroupB , 'Biological_Replicate', 'Different_Sample')

ggplot(dt.melt.cor, aes(x=Biological.Group.Type, y= cor))+
  geom_violin(trim = F, aes(fill = Biological.Group.Type), color = 'white')+
  geom_boxplot(width=.1, size = .5, color = 'gray40')+
  stat_compare_means()+
  labs(x= 'Group Type', 
       y = 'Spearman correlation')+
  theme_bw()+
  scale_fill_manual(values = c('#070d59', '#ceddef'))+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=20, face='bold',color='gray40'),
        axis.text.y = element_text(size = 20, face = 'bold',    color='gray40'),
        axis.title = element_text(size=20,face='bold'),
        legend.title = element_text(size=20,face='bold'),
        legend.text = element_text(size=20,face='bold',color='gray40'))
ggsave('./charts/Fig.S1c.Correlation_violin.pdf', width = 9, height = 9)

compare_means(cor~Biological.Group.Type, data = dt.melt.cor, method = "wilcox.test")


  



# statistic for plot ------------------------------------------------------

dt.melt.cor.rep <- dt.melt.cor[Biological.Group.Type=='Different_Sample']
median(dt.melt.cor.rep$cor)
sd(dt.melt.cor.rep$cor)

lst.cor.rep <- dt.melt.cor[Biological.Group.Type=='Biological_Replicate']$cor
lst.cor.diff <- dt.melt.cor[Biological.Group.Type=='Different_Sample']$cor
t.test <- t.test(lst.cor.diff, lst.cor.rep)
t.test$p.value
compare_means(cor~Biological.Group.Type, data = dt.melt.cor)

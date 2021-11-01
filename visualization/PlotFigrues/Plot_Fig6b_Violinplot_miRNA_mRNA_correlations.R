library(ggplot2)
library(magrittr)
library(data.table)
dt.cor.result <- readRDS('./RDS/Corr_miRNA_Gene_combined.rds')



# Plot total 6 types of correlations relationships ------------------------
# used in manuscript

# dt.cor.result.miRNA.gene <- dt.cor.result[Type%in%c('miRNA-all gene', 'miRNA-host gene (same)', 'miRNA-host gene (opp)')]

ggplot(dt.cor.result, aes(x=Type, y = corr))+
  geom_violin(trim = F, aes(fill=Type, color = Type))+
  geom_boxplot(outlier.alpha=0, width = 0.1)+
  scale_y_continuous(limits=c(-1,1))+
  scale_fill_manual(values=c('#D6A09E','#ADD4EF','#ADD4EF'))+
  scale_color_manual(values=c('#D6A09E','#ADD4EF','#ADD4EF'))+
  labs(y = 'Correlation coefficients', x = '')+
  theme_bw()+
  theme(axis.text.y=element_text(size=18, face='bold',color = 'black'),
        axis.text.x=element_text(size=18, face='bold',color = 'black'),
        axis.title.y=element_text(size=18,face='bold', color = 'black'),
        legend.title=element_text(size=18,face='bold', color = 'black'),
        legend.text=element_text(size=18,face='bold', color = 'black'),
        legend.position = 'top')

ggsave('./charts/Fig6b_Correlation_violin.pdf', width = 12, height = 6.5)



# stats -------------------------------------------------------------------
dt.cor.result$Type%>%unique()
mean(dt.cor.result[Type=="miRNA-all gene"]$corr)
sd(dt.cor.result[Type=="miRNA-all gene"]$corr)
mean(dt.cor.result[Type=="miRNA-host gene (same)"]$corr)
sd(dt.cor.result[Type=="miRNA-host gene (same)"]$corr)
mean(dt.cor.result[Type=="miRNA-host gene (oppo)"]$corr)
sd(dt.cor.result[Type=="miRNA-host gene (oppo)"]$corr)

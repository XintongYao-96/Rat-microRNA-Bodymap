library(data.table)
library(magrittr)
library(ggplot2)

dt.base.quality <- read.table('./tables/mqc_fastqc_per_base_sequence_quality_plot_1.txt')
setnames(dt.base.quality, c('Sample.ID', seq(1,50)))
dt.base.quality.forPlot <- dt.base.quality[-1, ]
dt.base.quality <- melt(dt.base.quality.forPlot, id.var = 'Sample.ID')


ggplot(dt.base.quality, aes(x = variable, y = value))+
  geom_boxplot(color = '#2A4318')+
  geom_hline(yintercept = 30, color = '#72B568', size = 2, linetype = 'dashed')+
  geom_hline(yintercept = 20, color = '#E8DBC3', size = 2, linetype = 'dashed')+
  scale_y_continuous(limits = c(0, 40))+
  labs(x=  'Base',
       y = 'Phred quality score')+
  scale_y_continuous(limits = c(25, 40))+
  theme_bw()+  
  theme(
    plot.title = element_text(size = 18, color = 'black', face = 'bold', hjust = 0.5),
    axis.text.x = element_text(size=12,color = 'black',face = 'bold', hjust = 0.5),
    axis.text.y = element_text(size = 14, color = 'black',face = 'bold'),
    axis.title = element_text(size=18,color = 'black',face='bold'),
    legend.title = element_text(size=18,color = 'black',face='bold'),
    legend.text = element_text(size=18,color = 'black',face='bold'))
ggsave('./charts/Fig2b_Boxplot_BaseQuality.pdf', width = 12, height = 6)  

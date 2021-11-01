library(data.table)
library(magrittr)
library(ggplot2)
source('./scripts/parameters.R')

dt.meta <- fread('./tables/Table_metadata_20210428.csv')
dt.forPlot <- dt.meta[, .(Sample_ID, total_sequences_before_trim)]
setnames(dt.forPlot, c('Sample_ID', 'Total_read'))
dt.forPlot$Total_reads_M <- dt.forPlot$Total_read/1000000


ggplot(dt.forPlot, aes(x= Sample_ID, y = Total_reads_M))+
  geom_bar(stat="identity", fill = color.miRNA.type[2],  width = 1)+
  scale_fill_manual(values = color.miRNA.type)+
  theme_classic()+  
  labs(y = 'Total number of sequencing reads (M)')+
  theme(
    legend.position = 'top',
    plot.title = element_text(size = 18, face = 'bold', hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18, face = 'bold'),
    axis.title = element_text(size=18,face='bold'),
    legend.title = element_text(size=18,face='bold'),
    legend.text = element_text(size=18,face='bold'))
ggsave('./charts/Fig2b_Barplot_total_sequence.pdf', width = 9, height = 5)

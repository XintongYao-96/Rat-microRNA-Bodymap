rm(list = ls())
library(data.table)
library(magrittr)
library(RColorBrewer)
library(ggplot2)

dt.readstats <- readRDS('./RDS_202201/ReadStats.rds')
dt.meta <- readRDS('./RDS/metadata.rds')

dt.readstats.forPlot <- dt.readstats[Type%in%c('Total_sequence', 'Total_forAlign')]
dt.readstats.forPlot[, Count_M := Count/1000000]
dt.readstats.dcast <- dcast(dt.readstats.forPlot, Sample.ID~Type, value.var  = 'Count_M')
dt.readstats.dcast[, Reads_filtered := Total_sequence-Total_forAlign]
dt.readstats.plot <- dt.readstats.dcast[, .(Sample.ID, Reads_filtered, Total_forAlign)]
setnames(dt.readstats.plot, c("Sample.ID",  "Reads_failed_filter", "Reads_passed_filter"))
dt.forPlot <- melt(dt.readstats.plot)%>%setnames(., c('Sample_ID', 'Type', 'Count'))


# color setting
color.type = brewer.pal(3, 'Set2')[c(2,1)]

ggplot(dt.forPlot, aes(x= Sample_ID, y = Count, fill = Type))+
  geom_bar(stat="identity",   width = 1)+
  scale_fill_manual(values = color.type)+
  theme_classic()+  
  labs(y = 'Count (M)')+
  theme(legend.position = 'top',
    plot.title = element_text(size = 18, face = 'bold', hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, face = 'bold', color = 'black'),
    axis.title = element_text(size = 14, face = 'bold', color = 'black'),
    legend.title = element_text(size = 14, face = 'bold', color = 'black'),
    legend.text = element_text(size = 14, face = 'bold', color = 'black'))
ggsave('./charts/Fig2a_Barplot_total_sequence.pdf', width = 9, height = 5)



# Total Genome ------------------------------------------------------------

dt.readstats.forPlot <- dt.readstats[Type%in%c('Total_forAlign', 'total_genome')]
dt.readstats.dcast <- dcast(dt.readstats.forPlot, Sample.ID~Type, value.var  = 'Count')
dt.readstats.dcast[, Genome_ratio := total_genome/Total_forAlign*100]

dt.readstats.dcast$Sample.ID <- factor(dt.readstats.dcast$Sample.ID, levels = dt.meta$Colnames.miRNA)


ggplot(dt.readstats.dcast, aes(x= Sample.ID, y = Genome_ratio))+
  geom_bar(stat="identity",   width = 1)+
  theme_classic()+  
  labs(y = 'Count (M)')+
  theme(legend.position = 'top',
        plot.title = element_text(size = 18, face = 'bold', hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, face = 'bold', color = 'black'),
        axis.title = element_text(size = 14, face = 'bold', color = 'black'),
        legend.title = element_text(size = 14, face = 'bold', color = 'black'),
        legend.text = element_text(size = 14, face = 'bold', color = 'black'))

rm(list = ls())
library(data.table)
library(magrittr)
library(RColorBrewer)
library(ggplot2)

dt.meta <- readRDS('./RDS/metadata.rds')
dt.readstats <- readRDS('./RDS_202201/ReadStats.rds')


# Sample filtering
IDs.sample.failedCluster <- c("lane3_Brn_M_6_3_1", "lane3_Spl_F_100_4_1")
# dt.readstats.passFilt <- dt.readstats[!Sample.ID%in%IDs.sample.failedCluster]


# Plot reads component ----------------------------------------------------

dt.readstats.comp <- dt.readstats[!Type%in%c('Total_sequence', 'Total_forAlign', 'total_genome')]
dt.readstats.comp[, Proportion := Count/sum(Count)*100, by = .(Sample.ID)]
# dt.readstats.comp[, check := sum(Perc), by = .(Sample.ID)]

# Readstas for Plot, annotated metadata, get the order correct

setkey(dt.readstats.comp, 'Sample.ID')
dt.readstats.comp.forPlot <- dt.readstats.comp[dt.meta$Colnames.miRNA] 
dt.readstats.comp.forPlot$Type <- factor(dt.readstats.comp.forPlot$Type, 
                                         levels = c("miRNA", "piRNA", "tRNA", "rRNA", "miscRNA","mRNA", "other_RNA","other_from_Genome", "unmapped"))
dt.readstats.comp.forPlot$Sample.ID <- factor(dt.readstats.comp.forPlot$Sample.ID,
                                              levels = dt.meta$Colnames.miRNA)
# Color setting

color.miRNA.type <- data.table(Type = c("miRNA", "piRNA", "tRNA", "rRNA", "miscRNA","mRNA", "other_RNA","other_from_Genome", "unmapped"),
                               Color = c("#9ECAE1", "#3182BD",'#9C9BC4', '#736BAE', '#4E288A', '#ACD8A1', '#5FA865', '#255A36','gray80')
)


ggplot(dt.readstats.comp.forPlot, aes(x= Sample.ID, y = Proportion, fill = Type))+
  geom_bar(stat="identity",position="stack",  width = 1)+
  scale_fill_manual(values = color.miRNA.type$Color)+
  theme_classic()+  
  labs(y = 'Proportion (%)')+
  theme(
    plot.title = element_text(size = 18, face = 'bold', hjust = 0.5),
    axis.text.x = element_text(size = 4, face = 'bold', color = 'black', angle = 90),
    axis.text.y = element_text(size = 18, face = 'bold',color = 'black'),
    axis.title = element_text(size=18,face='bold'),
    legend.title = element_text(size=18,face='bold'),
    legend.text = element_text(size=18,face='bold'))
ggsave('./charts/Fig2b_Barplot_ReadsComponent.pdf', width = 16, height = 8)



# stats -------------------------------------------------------------------

mean(dt.readstats[Type=='Total_forAlign']$Count)
sd(dt.readstats[Type=='Total_forAlign']$Count)

IDs.mature.tst <- dt.meta[Organ=='Tst'][Age%in%c(6, 21)]$Colnames.miRNA

dt.readstats.comp.forPlot[Sample.ID%in%IDs.mature.tst][Type=='piRNA']$Proportion%>%mean()
dt.readstats.comp.forPlot[Sample.ID%in%IDs.mature.tst][Type=='piRNA']$Proportion%>%sd()

dt.readstats.comp.forPlot[!Sample.ID%in%IDs.mature.tst][Type=='piRNA']$Proportion%>%mean()
dt.readstats.comp.forPlot[!Sample.ID%in%IDs.mature.tst][Type=='piRNA']$Proportion%>%sd()


dt.readstats.comp.forPlot[!Sample.ID%in%IDs.mature.tst][Type=='miRNA']$Proportion%>%mean()
dt.readstats.comp.forPlot[!Sample.ID%in%IDs.mature.tst][Type=='miRNA']$Proportion%>%sd()

dt.readstats.comp.forPlot[!Sample.ID%in%IDs.mature.tst][Type=='piRNA']$Proportion%>%mean()
dt.readstats.comp.forPlot[!Sample.ID%in%IDs.mature.tst][Type=='piRNA']$Proportion%>%sd()

dt.map <-dt.readstats.comp.forPlot[Type=='unmapped']
dt.map[, mapped := 100- Proportion]
dt.map$mapped
mean(dt.map$Proportion)
sd(dt.map$mapped)

0.0098*1000000
dt.readstats.comp.forPlot[Type=='miRNA'][Proportion<25][]

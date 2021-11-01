library(data.table)
library(magrittr)
library(ggplot2)
source('./scripts/PlotFigures/parameters.R')

dt.meta <- readRDS('./RDS/metadata_with_miRNA_Number.rds')

dt.meta.miRNA.number <- dt.meta[, .(ID.sample, anot_miRNA_number, novel_miRNA_number)]
setnames(dt.meta.miRNA.number, c('Sample_ID', 'annotated miRNA', 'novel miRNA'))

dtforPlot <- melt(dt.meta.miRNA.number, id.vars = 'Sample_ID',
                   variable.name = 'miRNA.type',
                   value.name = 'Number')
dtforPlot$miRNA.type <- factor(dtforPlot$miRNA.type, levels = rev(c(  'annotated miRNA', 'novel miRNA') ))

color.miRNA.type <- c('#AED3EF', '#2D679B')

ggplot(dtforPlot, aes(x= Sample_ID, y = Number, fill = miRNA.type))+
  geom_bar(stat="identity",position="stack",  width = 1)+
  scale_fill_manual(values = color.miRNA.type)+
  theme_classic()+  
  labs(y = 'Number of miRNAs detected')+
  theme(
    legend.position = 'top',
    plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
    axis.text.x = element_text(size = 14, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 14, face = 'bold',color = 'black'),
    axis.title = element_text(size=14,face='bold'),
    legend.title = element_text(size=14,face='bold'),
    legend.text = element_text(size=14,face='bold'))
ggsave('./charts/Fig7a_Barplot_miRNA_number.pdf', width = 9, height = 5)  


# stats ------------------------------------------------------------------

dt.meta <- readRDS('./RDS/metadata_with_miRNA_Number.rds')
mean(dt.meta[Organ=='Brn']$anot_miRNA_number)
sd(dt.meta[Organ=='Brn']$anot_miRNA_number)

mean(dt.meta[Organ=='Brn']$novel_miRNA_number)
sd(dt.meta[Organ=='Brn']$novel_miRNA_number)

mean(dt.meta[Organ=='Lvr']$anot_miRNA_number)
sd(dt.meta[Organ=='Lvr']$anot_miRNA_number)

mean(dt.meta[Organ=='Lvr']$novel_miRNA_number)
sd(dt.meta[Organ=='Lvr']$novel_miRNA_number)


dt.meta[, .(N.novel.mean = mean(novel_miRNA_number)), by = .(Organ)]
dt.meta[, .(N.novel.sd = sd(novel_miRNA_number)), by = .(Organ)]
dt.meta[, .(N.known.mean = mean(anot_miRNA_number)), by = .(Organ)]
dt.meta[, .(N.known.sd = sd(anot_miRNA_number)), by = .(Organ)]

# 
# library(pheatmap)
# 
# df.annot.sample.id <- data.frame(Sample.ID = dt.meta$Sample_ID,
#                                  sex = dt.meta$Sex, row.names = 1,
#                                  Age = dt.meta$Age,
#                                  Organ = dt.meta$Organ_Abbr_3chars)
# 
# 
# annot_color <- list(Organ=colors.organ,
#                     Age=colors.age,
#                     Sex=color.sex )
# 
# 
# mat <- data.frame(sample.id =  dt.meta$Sample_ID,
#                   r1 = rep(1,320),
#                   r2 = rep(2,320), row.names = 1)%>%as.matrix()%>%t()
# 
# pheatmap(mat,
#          cluster_rows = F, cluster_cols = F,
#          annotation_col = df.annot.sample.id,
#          annotation_colors = annot.color,
#          filename = './charts/Fig3b_annot_bar.pdf', width = 14, height = 8)

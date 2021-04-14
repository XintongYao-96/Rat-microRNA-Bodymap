library(data.table)
library(magrittr)
library(ggplot2)
library(Rtsne)
library(ggpubr)
library(ggthemes)
source('scripts/parameters.R')
dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
dt.meta <- readRDS('./RDS/metadata.rds')[ID.sample%in%colnames(dt.exprmat)]

# tSNE
set.seed(1)
tsne.info <- Rtsne(t(dt.exprmat), perplexity = 32)
dt.tSNE <- data.table(tsne.info$Y)
colnames(dt.tSNE) <- c('tSNE_1', 'tSNE_2')
dt.tSNE$ID.sample<- colnames(dt.exprmat)


# dtforPLot 

dt.forPlot <- merge(dt.tSNE,dt.meta, by= 'ID.sample')

ggscatter(dt.forPlot, x='tSNE_1', 'tSNE_2', color = 'Organ', shape = 'Age',
          ellipse = F,
          size = 3)+
  scale_color_manual(values = colors.organ)+
  theme_base()
ggsave('./charts/Fig.S1c.tSNE.pdf', width = 8, height = 7)

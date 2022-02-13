rm(list = ls())

library(data.table)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
source('./scripts/PlotFigures/parameters.R')


dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.forPCA <- readRDS('./RDS/exprMat/exprMat_logCPM_r604c318.rds')
dt.meta.forPCA <- dt.meta[ID.sample%in%colnames(exprMat.forPCA)]


# Perform PCA -------------------------------------------------------------

pca_prcomp <- prcomp(t(exprMat.forPCA), scale. = F)
pcs <- predict(pca_prcomp)%>%as.data.frame()
pcs$ID.sample <- rownames(pcs)
dim(pcs)
Mat.PC.loadings <- pca_prcomp$rotation

dt.meta.pcs <- merge(dt.meta, pcs, by ='ID.sample', by.all=T)%>%as.data.table()

saveRDS(dt.meta.pcs, './RDS/meta_PCS.rds')
saveRDS(Mat.PC.loadings, './RDS/PC_loadings.rds')


# set colors --------------------------------------------------------------
# 
# lst.organ.age <- unique(dt.meta.pcs$organ_age )
# 
# lst.color.Adr <- brewer.pal(9, 'Blues')[c(3,5,7,9)]
# lst.color.Spl <- brewer.pal(9, 'BuGn')[c(3,5,7,9)]
# lst.color.Hrt <- brewer.pal(9, 'Reds')[c(3,5,7,9)]
# lst.color.Kdn <- brewer.pal(9, 'YlOrRd')[c(3,5,7,9)]
# lst.color.Lng <- brewer.pal(9, 'YlOrBr')[c(3,5,7,9)]
# lst.color.Lvr <- brewer.pal(9, 'RdPu')[c(3,5,7,9)]
# lst.color.Msc <- brewer.pal(9, 'Purples')[c(3,5,7,9)]
# lst.color.Brn <- brewer.pal(9, 'PuRd')[c(3,5,7,9)]
# lst.color.Thm <- brewer.pal(9, 'PuBuGn')[c(3,5,7,9)]
# lst.color.Tst <- brewer.pal(9, 'Oranges')[c(3,5,7,9)]
# lst.color.Utr <- brewer.pal(9, 'BuPu')[c(3,5,7,9)]
# 
# color.organ.age <- c(lst.color.Adr,lst.color.Brn,lst.color.Hrt,lst.color.Kdn,lst.color.Lng,lst.color.Lvr,lst.color.Msc,lst.color.Spl,
#                      lst.color.Thm,lst.color.Tst,lst.color.Utr)
# names(color.organ.age) <- lst.organ.age



# Plot PCA ----------------------------------------------------------------

#PC1-PC2
ggplot(dt.meta.pcs, aes(x=PC1, y=PC2))+
  geom_point(aes(color=Organ, shape = Age),size=2.5)+
  scale_color_manual(values = colors.organ)+
  theme_bw()+
  labs(x=sprintf("PC1 (%.1f%%)", summary(pca_prcomp)$importance[2,1]*100),
       y=sprintf("PC2 (%.1f%%)", summary(pca_prcomp)$importance[2,2]*100,1),
       color='Sample')+
  theme(legend.position = 'none',
        axis.text = element_text(size=16,color='black'),
        axis.title = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16,color='black'))
ggsave('./charts/Fig4b_PCA_PC12.pdf', width = 6, height = 6)

# pc2-pc3
ggplot(dt.meta.pcs, aes(x=PC2, y=PC3))+
  geom_point(aes(color=Organ, shape = Age),size=2.5)+
  scale_color_manual(values = colors.organ)+
  theme_bw()+
  labs(x=sprintf("PC2 (%.1f%%)", summary(pca_prcomp)$importance[2,2]*100),
       y=sprintf("PC3 (%.1f%%)", summary(pca_prcomp)$importance[2,3]*100,1),
       color='Sample')+
  theme(legend.position = 'none',
         axis.text = element_text(size=16,color='black'),
         axis.title = element_text(size=16),
         legend.title = element_text(size=16),
         legend.text = element_text(size=16,color='black'))
ggsave('./charts/Fig4c.PCA_PC23.pdf', width = 6, height =6)

#PC1-PC3
  ggplot(dt.meta.pcs, aes(x=PC1, y=PC3))+
    geom_point(aes(color=Organ, shape = Age),size=2.5)+
    scale_color_manual(values = colors.organ)+
    theme_bw()+
    labs(x=sprintf("PC1 (%.1f%%)", summary(pca_prcomp)$importance[2,1]*100),
         y=sprintf("PC3 (%.1f%%)", summary(pca_prcomp)$importance[2,2]*100,1),
         color='Sample')+
    theme(legend.position = 'none',
          axis.text = element_text(size=16,color='black'),
          axis.title = element_text(size=16),
          legend.title = element_text(size=16),
          legend.text = element_text(size=16,color='black'))
ggsave('./charts/Fig4d.PCA_PC13.pdf', width = 6, height = 6)
  

# test --------------------------------------------------------------------

dt.meta.pcs$size = 1
dt.meta.pcs[Colnames.miRNA=='lane5_Brn_F_100_4_1']$size = 2
ggplot(dt.meta.pcs, aes(x=PC1, y=PC3))+
  geom_point(aes(color=Organ, shape = Age, size = size))+
  scale_color_manual(values = colors.organ)+
  theme_bw()+
  labs(x=sprintf("PC1 (%.1f%%)", summary(pca_prcomp)$importance[2,1]*100),
       y=sprintf("PC3 (%.1f%%)", summary(pca_prcomp)$importance[2,2]*100,1),
       color='Sample')+
  theme(legend.position = 'none',
        axis.text = element_text(size=16,color='black'),
        axis.title = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16,color='black'))

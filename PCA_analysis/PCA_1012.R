source('./scripts/library.R')
source('./scripts/parameters.R')


dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.forPCA <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
dt.meta.forPCA <- dt.meta[ID.sample%in%colnames(exprMat.forPCA)]

pca_prcomp <- prcomp(t(exprMat.forPCA), scale. = F)
pcs <- predict(pca_prcomp)%>%as.data.frame()
pcs$ID.sample <- rownames(pcs)
dim(pcs)
Mat.PC.loadings <- pca_prcomp$rotation

dt.meta.pcs <- merge(dt.meta, pcs, by='ID.sample', by.all=T)%>%as.data.table()

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

# pc2-pc3
ggplot(dt.meta.pcs, aes(x=PC2, y=PC3))+
  geom_point(aes(color=Organ, shape = Age),size=5)+
  scale_color_manual(values = colors.organ)+
  theme_bw()+
  labs(x=sprintf("PC2 (%.1f%%)", summary(pca_prcomp)$importance[2,2]*100),
       y=sprintf("PC3 (%.1f%%)", summary(pca_prcomp)$importance[2,3]*100,1),
       color='Sample')+
  theme(axis.text = element_text(size=12,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
  ggsave('./charts/Fig1b.PCA.png', width = 9, height = 8)

#PC1-PC2
ggplot(dt.meta.pcs, aes(x=PC1, y=PC2))+
  geom_point(aes(color=Organ, shape = Age),size=5)+
  scale_color_manual(values = colors.organ)+
  theme_bw()+
  labs(x=sprintf("PC1 (%.1f%%)", summary(pca_prcomp)$importance[2,1]*100),
       y=sprintf("PC2 (%.1f%%)", summary(pca_prcomp)$importance[2,2]*100,1),
       color='Sample')+
  theme(axis.text = element_text(size=12,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
  ggsave('./charts/FigS1.PCA1_2.png', width = 9, height = 8)

#PC1-PC3
ggplot(dt.meta.pcs, aes(x=PC1, y=PC3))+
  geom_point(aes(color=organ_age),size=3)+
  scale_color_manual(values = color.organ.age)+
  theme_bw()+
  labs(x=sprintf("PC1 (%.1f%%)", summary(pca_prcomp)$importance[2,1]*100),
       y=sprintf("PC3 (%.1f%%)", summary(pca_prcomp)$importance[2,3]*100,1),
       color='Sample')+
  theme(axis.text = element_text(size=12,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))+
  ggsave('./charts/Plot_PCA_pc1_3.png', width = 14, height = 10)

#PC3-PC4
ggplot(dt.meta.pcs, aes(x=PC3, y=PC4))+
  geom_point(aes(color=organ_age),size=3)+
  scale_color_manual(values = color.organ.age)+
  theme_bw()+
  labs(x=sprintf("PC3 (%.1f%%)", summary(pca_prcomp)$importance[2,3]*100),
       y=sprintf("PC4 (%.1f%%)", summary(pca_prcomp)$importance[2,4]*100,1),
       color='Sample')+
  theme(axis.text = element_text(size=12,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))+
  ggsave('./charts/Plot_PCA_pc3_4.png', width = 14, height = 10)


ggplot(dt.meta.pcs, aes(x=PC3, y=PC4))+
  geom_point(aes(color=organ_age),size=3)+
  scale_color_manual(values = color.organ.age)+
  theme_bw()+
  labs(x=sprintf("PC3 (%.1f%%)", summary(pca_prcomp)$importance[2,3]*100),
       y=sprintf("PC4 (%.1f%%)", summary(pca_prcomp)$importance[2,4]*100,1),
       color='Sample')+
  theme(axis.text = element_text(size=12,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
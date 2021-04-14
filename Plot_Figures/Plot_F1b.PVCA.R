rm(list = ls())
library(data.table)
source("./scripts/Functions/Function.Calc.PVCA.R")



dt.meta <- readRDS('./RDS/metadata.rds')
setkey(dt.meta,ID.sample)



# PVCA analysis of miRNA profile
exprMat.miRNA <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
pvca.miRNA <- Function.calc.PVCA(exprMat=exprMat.miRNA,
                                 expDesign = dt.meta[,.(ID.sample,Organ,Age,Sex)])
write.csv(pvca.miRNA,'./tables/PVCA of miRNA.csv',quote=FALSE,row.names=FALSE)

rm(exprMat.miRNA)
gc()


# PVCA analysis of mRNA profile
exprMat.mRNA <- readRDS('./RDS/exprMat/exprMat_mRNA_log2FPKM.rds')
pvca.mRNA <-  Function.calc.PVCA(exprMat=exprMat.mRNA,
                                 expDesign = dt.meta[,.(ID.sample,Organ,Age,Sex)])
write.csv(pvca.mRNA,'./tables/PVCA of mRNA.csv',quote=FALSE,row.names=FALSE)

rm(exprMat.mRNA)
gc()


# Bind together
pvca.miRNA$Type='miRNA'
pvca.mRNA$Type='gene'
pvca.combine <- rbind(pvca.miRNA,pvca.mRNA)


library(ggplot2)
# Old version
# ggplot(pvca.combine,aes(y=reorder(Stage, Perc.effect),x=Perc.effect))+
#   geom_bar(aes(fill=Type),stat = 'identity', pos='dodge',width=.8)+
#   labs(y='',x='Percentage of overall variance')+
#   theme_bw()+
#   theme(axis.text.y=element_text(size=18,face='bold'),
#         axis.text.x=element_text(size=18,hjust=1,face='bold'),
#         axis.title.y=element_text(size=20,face='bold'),
#         legend.title=element_text(size=20,face='bold'),
#         legend.text=element_text(size=18,face='bold'))+
#   scale_fill_brewer(palette = 'Set1')
# 
# ggsave('./charts/Fig.1b.PVCA.pdf',width=6,height=7)
# 

ggplot(pvca.combine,aes(x=reorder(Stage, -Perc.effect),y=Perc.effect))+
  geom_bar(aes(fill=Type),stat = 'identity', pos='dodge',width=.8)+
  labs(x='',y='Percentage of overall variance')+
  theme_bw()+
  theme(axis.text.y=element_text(size=18,face='bold'),
        axis.text.x=element_text(size=18,face='bold'),
        axis.title.y=element_text(size=20,face='bold'),
        legend.title=element_text(size=20,face='bold'),
        legend.text=element_text(size=18,face='bold'))+
  scale_fill_manual(values = c('#2D679B','#D36A87'))
ggsave('./charts/Fig.1b.PVCA.pdf',width=12,height=6.5)




pvca.combine.dcast <- dcast(pvca.combine, Stage~Type, value.var = 'Perc.effect')

rm(list = ls())
library(data.table)
library(ggplot2)
source("./scripts/Functions/Function.Calc.PVCA.R")
dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.miRNA <- readRDS('./RDS/exprMat/exprMat_logCPM_r604c318.rds')
exprMat.mRNA <- readRDS('./RDS/exprMat/exprMat_mRNA_log2FPKM.rds')




# PVCA --------------------------------------------------------------------


# PVCA analysis of miRNA profile
setkey(dt.meta,ID.sample)

pvca.miRNA <- Function.calc.PVCA(exprMat=exprMat.miRNA,
                                 expDesign = dt.meta[,.(ID.sample,Organ,Age,Sex)])
# write.csv(pvca.miRNA,'./tables/PVCA of miRNA.csv',quote=FALSE,row.names=FALSE)
rm(exprMat.miRNA)
gc()


# PVCA analysis of mRNA profile
pvca.mRNA <-  Function.calc.PVCA(exprMat=exprMat.mRNA,
                                 expDesign = dt.meta[,.(ID.sample,Organ,Age,Sex)])
# write.csv(pvca.mRNA,'./tables/PVCA of mRNA.csv',quote=FALSE,row.names=FALSE)
rm(exprMat.mRNA)
gc()


# Bind together
pvca.miRNA$Type='miRNA'
pvca.mRNA$Type='gene'
pvca.combine <- rbind(pvca.miRNA,pvca.mRNA)



# plot --------------------------------------------------------------------


ggplot(pvca.combine,aes(x=reorder(Stage, -Perc.effect),y=Perc.effect))+
  geom_bar(aes(fill=Type),stat = 'identity', pos='dodge',width=.8)+
  labs(x='',y='Percentage of overall variance (%)')+
  theme_bw()+
  theme(axis.text.y=element_text(size=18, face='bold',color = 'black'),
        axis.text.x=element_text(size=18, face='bold',color = 'black'),
        axis.title.y=element_text(size=18,face='bold', color = 'black'),
        legend.title=element_text(size=18,face='bold', color = 'black'),
        legend.text=element_text(size=18,face='bold', color = 'black'),
        legend.position = 'top')+
  scale_fill_manual(values = c('#2D679B','#D36A87'))
ggsave('./charts/Fig6a_PVCA.pdf',width=12,height=6.5)




pvca.combine.dcast <- dcast(pvca.combine, Stage~Type, value.var = 'Perc.effect')


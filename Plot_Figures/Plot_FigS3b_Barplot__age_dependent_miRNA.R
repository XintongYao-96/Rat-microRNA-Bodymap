rm(list = ls())
library(data.table)
library(magrittr)

source('./scripts/parameters.R')
dt.DEGs.forw <- readRDS('./RDS/DEG_results/DEG_AgeperOrgan.rds')


dt.DEGs.rev <- data.table(
  logFC = -dt.DEGs.forw$logFC,
  AveExpr = dt.DEGs.forw$AveExpr,
  t = dt.DEGs.forw$t,
  P.Value = dt.DEGs.forw$P.Value,
  adj.P.Val = dt.DEGs.forw$adj.P.Val,
  B = dt.DEGs.forw$B,
  ID.gene = dt.DEGs.forw$ID.gene,
  Age.A = dt.DEGs.forw$Age.B,
  Age.B = dt.DEGs.forw$Age.A,
  Organ = dt.DEGs.forw$Organ
)

dt.DEGs <- rbindlist(list(dt.DEGs.forw,dt.DEGs.rev))

rm(dt.DEGs.forw,dt.DEGs.rev)


dt.DEGs[,Group:=sprintf("%svs%s",Age.A,Age.B)]


dt.DEGs[,adj.P.Val:=P.Value*6]
dt.DEGs[,Type.regulate:=ifelse(adj.P.Val<0.05,ifelse(logFC>0,'Up','Down'),'Maintain')]


lst.develop.depend.miRNA <- dt.DEGs[Type.regulate!='Maintain']$ID.gene%>%unique()
length(lst.develop.depend.miRNA)
saveRDS(lst.develop.depend.miRNA, './RDS/list_development_dependent_miRNA.rds')

# Fig. 4a | Barplot of number of DEGs in each development stage -----------


dt.DEGs.forPlot <- dt.DEGs[Group%in%c('2vs6','2vs21','2vs104')][adj.P.Val<0.05]

ct.DEGs <- dt.DEGs.forPlot[,.N,by=.(Organ,Group,Type.regulate)]
ct.DEGs[,Y:=ifelse(Type.regulate=='Up',N,-N)]


ct.DEGs[N==min(ct.DEGs$N)]
ct.DEGs[N==max(ct.DEGs$N)]

library(ggplot2)
library(RColorBrewer)

ct.DEGs$Type.regulate <- factor(ct.DEGs$Type.regulate,levels=names(colors.regulateType))
ct.DEGs$Group <- factor(ct.DEGs$Group,levels=c('2vs6','2vs21','2vs104'))

ggplot(ct.DEGs,aes(x=Group,y=Y))+
  geom_bar(stat = 'identity',pos='stack',color='gray80',
           aes(fill=Type.regulate))+
  # geom_point()+
  theme_bw()+
  labs(x='',y='Number of development-dependent miRNAs',fill='')+
  theme(axis.text.x=element_text(size=12,face='bold',angle=50,hjust=1),
        axis.text.y=element_text(size=10,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.position = c(0.92,0.2),
        legend.text = element_text(size=12,face='bold'),
        strip.background = element_rect(color='gray20',fill='gray95'),
        strip.text = element_text(size=14,face='bold')
  )+
  scale_fill_manual(values=brewer.pal(3,'Set2')[c(2,1)])+
  facet_wrap(~Organ,nrow = 2)

ggsave('./charts/Fig.S3b.Numbers_of_DEG_in_each_development_per_organ.pdf',width=10,height=6)  

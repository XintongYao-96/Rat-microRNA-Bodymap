rm(list = ls())
source('./scripts/library.R')
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


# Figure S3| Barplot of number of DEGs in each development stage -----------


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
  labs(x='',y='Number of develop-dependent miRNAs',fill='')+
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

ggsave('./charts/Fig.4a.Numbers of DEGs in each development per organ.png',width=10,height=5)  






# Figure 4b organ pattern proportion---------------------------------------------------------------

dt.DEGs.forPlot <- dt.DEGs[Group%in%c('2vs6','6vs21','21vs104')]

abbrs.Type.regulate <- c(Up='U',Down='D',Maintain='M')
dt.DEGs.forPlot[,Abbr.Type.Regulate:=abbrs.Type.regulate[Type.regulate]]

dt.pattern <- dcast(dt.DEGs.forPlot, ID.gene+Organ~Group, value.var='Abbr.Type.Regulate',fill='M')
dt.pattern[,Pattern:=paste(`2vs6`,`6vs21`,`21vs104`,sep='')]

ct.PatternPerOrgan <- dt.pattern[,.N,by=.(Pattern,Organ)]
ct.PatternPerOrgan[,Perc:=N/sum(N)*100,by=Organ]

ct.Pattern.dcast <- dcast(ct.PatternPerOrgan,Organ~Pattern,value.var='Perc',fill=0.00)


# add mean of each pattern to the plot
result <- data.frame( t(c('mean', apply(ct.Pattern.dcast[, 2:28], 2, mean)%>%as.numeric())))
names(result) <- names(ct.Pattern.dcast)
ct.Pattern.dcast <- data.table(rbind(ct.Pattern.dcast, result))
  


#write.csv(ct.Pattern.dcast,'./results/final/tables/Development pattern of miRNA in each organ.csv',quote=FALSE,row.names=FALSE)

dt.forPlot <- melt(ct.Pattern.dcast,id.vars = 'Organ',variable.name = 'Pattern', value.name = 'Perc')

dt.forPlot$Organ <- factor(dt.forPlot$Organ, levels = c('mean', rev(names(colors.organ))))
dt.forPlot$Perc <-  as.numeric(dt.forPlot$Perc)


ggplot(dt.forPlot,aes(x=Pattern,y=Organ))+
  geom_tile(aes(fill=Perc^(1/4)),color='gray20')+
  geom_text(aes(label=sprintf("%.2f",Perc)),fontface='bold')+
  labs(x='',y='')+
  theme_bw()+
  theme(axis.text.y=element_text(size=14,face='bold'),
        axis.text.x=element_text(size=10,face='bold'),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),)+
  scale_fill_gradient2(low ='#93abd3', mid = 'white', high = '#ec5858', midpoint = 0.8 )+
  scale_size_continuous(range=c(2,4))+
  scale_color_manual(values=c('white','gray80'))+
  guides(size=FALSE,color=FALSE,fill=FALSE)


ggsave('./charts/Fig_4b_PatternOfmiRNAInDevelopment.pdf',width=15,height=5)
saveRDS(dt.forPlot[Organ!='mean'], './RDS/miRNA_development_pattern_proportion.rds')



# explore trajectories of novel miRNAs ------------------------------------

dt.novel.miRNA <- readRDS('./RDS/novel_miRNA_analysis/miRNA_info.rds')
lst.novel.miRNA <- dt.novel.miRNA[type.miRNA=='novel']$miRNA.ID
dt.pattern.novel <- dt.pattern[ID.gene%in%lst.novel.miRNA]


ct.PatternPerOrgan <- dt.pattern.novel[,.N,by=.(Pattern,Organ)]
ct.PatternPerOrgan[,Perc:=N/sum(N)*100,by=Organ]

ct.Pattern.dcast <- dcast(ct.PatternPerOrgan,Organ~Pattern,value.var='Perc',fill=0.00)

dt.forPlot <- melt(ct.Pattern.dcast,id.vars = 'Organ',variable.name = 'Pattern', value.name = 'Perc')

dt.forPlot$Organ <- factor(dt.forPlot$Organ, levels = rev(names(colors.organ)))
dt.forPlot$Perc <-  as.numeric(dt.forPlot$Perc)

ggplot(dt.forPlot,aes(x=Pattern,y=Organ))+
  geom_tile(aes(fill=Perc^(1/4)),color='gray20')+
  geom_text(aes(label=sprintf("%.2f",Perc)),fontface='bold')+
  labs(x='',y='')+
  theme_bw()+
  theme(axis.text.y=element_text(size=14,face='bold'),
        axis.text.x=element_text(size=10,face='bold'),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),)+
  scale_fill_gradient2(low ='#93abd3', mid = 'white', high = '#ec5858', midpoint = 0.8 )+
  scale_size_continuous(range=c(2,4))+
  scale_color_manual(values=c('white','gray80'))+
  guides(size=FALSE,color=FALSE,fill=FALSE)
ggsave('./charts/Novel_PatternOfmiRNAInDevelopment.pdf',width=15,height=5)





# Print head (schematic diagram)

# 
# Stage2Int <- function(x){
#   if(x=='U') return (1)
#   if(x=='M') return (0)
#   if(x=='D') return (-1)
# }
# 
# dt.pattern.uni <- data.table(Pattern=colnames(ct.Pattern.dcast)[-1],S0=0)
# dt.pattern.uni$S1 <- dt.pattern.uni$S0 + sapply(dt.pattern.uni$Pattern, function(x)Stage2Int(substr(x,1,1)))
# dt.pattern.uni$S2 <- dt.pattern.uni$S1 + sapply(dt.pattern.uni$Pattern, function(x)Stage2Int(substr(x,2,2)))
# dt.pattern.uni$S3 <- dt.pattern.uni$S2 + sapply(dt.pattern.uni$Pattern, function(x)Stage2Int(substr(x,3,3)))
# dt.pattern.uni.melt <- melt(dt.pattern.uni,id.vars = 'Pattern',variable.name = 'Stage',value.name = 'y')
# dt.pattern.uni.melt[,y.scale:=(y-min(y))/(max(y)-min(y)),by=Pattern]
# dt.pattern.uni.melt$y.scale <- sapply(dt.pattern.uni.melt$y.scale,function(x)ifelse(is.na(x),1,x))
# dt.pattern.uni.melt[,x:=as.numeric(Stage)-1]
# 
# ggplot(dt.pattern.uni.melt,aes(x=x,y=y.scale,group=Pattern))+
#   geom_line(size=3,color='gray40',alpha=.6)+
#   geom_point(size=5,color='gray40',alpha=.6)+
#   labs(x='',y='')+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(),
#         panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())+
#   scale_y_continuous(limits=c(-0.5,1.5))+
#   facet_wrap(~Pattern,nrow = 1)
# 
# ggsave('./charts/Fig.4bHead.png',width=27,height=2)                  

# Print Head （profile)

patterns <- colnames(ct.Pattern.dcast)[-1]

exprMat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')

dt.exprSet <- as.data.table(melt(exprMat))
setnames(dt.exprSet,1:3,c('ID.gene','ID.sample','log2.CPM'))

dt.exprSet.Annot <- merge(dt.exprSet,dt.meta,by='ID.sample')
avg.exprSet <- dt.exprSet.Annot[,.(Avg=mean(log2.CPM)),by=.(ID.gene,Age,Organ)]
setkeyv(avg.exprSet,c('ID.gene','Organ'))

currPat<-patterns[3] 
pairs.miR.Organ <- dt.pattern[Pattern==currPat,.(ID.gene,Organ)]

avg.exprSet.forPlot <- merge(avg.exprSet,dt.pattern.novel,by=c('ID.gene','Organ'))
avg.exprSet.forPlot[,Avg.scale:=(Avg-mean(Avg))/sd(Avg),by=.(ID.gene,Organ,Pattern)]

ggplot(avg.exprSet.forPlot,aes(x=Age,y=Avg.scale,group=paste(Organ,ID.gene)))+
  geom_line(color='gray60',alpha=.1)+
  # geom_text(aes(label=Pattern),y=0)+
  labs(x='',y='')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face='bold',size=10),
        # axis.text = element_blank(),
        # axis.ticks = element_blank()
        )+
  facet_wrap(~Pattern, nrow=3)
  # facet_wrap(~Pattern,nrow = 1)

ggsave('./charts/Fig.4b.Head2_3_9.pdf',width=13,height=5)   

#########################
# pattern of novel miRNA#
#   20210206            #
#########################


avg.exprSet.forPlot <- merge(avg.exprSet,dt.pattern.novel,by=c('ID.gene','Organ'))
avg.exprSet.forPlot[,Avg.scale:=(Avg-mean(Avg))/sd(Avg),by=.(ID.gene,Organ,Pattern)]

ggplot(avg.exprSet.forPlot,aes(x=Age,y=Avg.scale,group=paste(Organ,ID.gene)))+
  geom_line(color='gray60',alpha=.1)+
  # geom_text(aes(label=Pattern),y=0)+
  labs(x='',y='')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face='bold',size=10),
        # axis.text = element_blank(),
        # axis.ticks = element_blank()
  )+
  facet_wrap(~Pattern, nrow=3)
ggsave('./charts/Novel_Head2_3_9.pdf',width=13,height=5)   


# facet_wrap(~Pattern,nrow = 1)


# generate numbers for manuscript -----------------------------------------
table(avg.exprSet.forPlot$Pattern)
8816/sum(table(avg.exprSet.forPlot$Pattern))
20/sum(table(avg.exprSet.forPlot$Pattern))
12/sum(table(avg.exprSet.forPlot$Pattern))
number.of.development.dependent <- length(unique(dt.DEGs.forPlot$ID.gene))
summary(dt.DEGs.forPlot$ID.gene)
ct.Pattern.dcast$MMM%>%mean()
ct.Pattern.dcast[Organ!='Thm']$MUU%>%mean()
ct.Pattern.dcast$UDU%>%mean()
ct.Pattern.dcast$DUD%>%mean()
ct.Pattern.dcast$DUU%>%mean()
ct.Pattern.dcast$DUU%>%mean()

avg.exprSet.forPlot[Pattern=='DDD']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='DDM']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='DDU']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='DMD']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='DMM']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='DMU']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='DUD']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='DUM']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='DUU']$ID.gene%>%unique()%>%length()


avg.exprSet.forPlot[Pattern=='MDD']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='MDM']$ID.gene%>%unique()%>%length()
avg.exprSet.forPlot[Pattern=='MDU']$ID.gene%>%unique()%>%length()


lst.pattern <- unique(avg.exprSet.forPlot$Pattern)
dt.pattern.n.miRNA <- sapply(lst.pattern, function(x){
                              pattern <- x
                              number.miRNA <-avg.exprSet.forPlot[Pattern== pattern]$ID.gene%>%unique()%>%length()
                              
                              dt.pattern.number <- data.table(pattern=pattern,
                                                              number.miRNA = number.miRNA)
                              
                              return(list(dt.pattern.number))
                            })%>%rbindlist()


# pathway enrichment analysis of UUU and DDD 1223------------------------------


lst.miRNA.uuu <- avg.exprSet.forPlot[Pattern=='UUU']$ID.gene%>%unique()
lst.miRNA.ddd <- avg.exprSet.forPlot[Pattern=='DDD']$ID.gene%>%unique()

dt.miRTarGene <- readRDS('./RDS/dt_MTI_update.rds')

# get target gene for cluster1&2miRNA

dt.uuu.miRTarGene <- dt.miRTarGene[ID.miRNA%in%lst.miRNA.uuu]
dt.ddd.miRTarGene <- dt.miRTarGene[ID.miRNA%in%lst.miRNA.ddd]
lst.uuu.miRTarGene <- dt.uuu.miRTarGene$ID.gene%>%unique()
lst.ddd.miRTarGene <- dt.ddd.miRTarGene$ID.gene%>%unique()

library(clusterProfiler)
library(org.Rn.eg.db)

Pathway.uuu.targene <- data.frame(enrichGO(lst.uuu.miRTarGene, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))
Pathway.ddd.gene <- data.frame(enrichGO(lst.ddd.miRTarGene, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))


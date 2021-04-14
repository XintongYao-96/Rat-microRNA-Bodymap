library(data.table)
library(magrittr)
library(ggplot2)
library(ggrepel)
source('./scripts/parameters.R')
dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
dt.meta <- readRDS('./RDS/metadata.rds')[ID.sample%in%colnames(dt.exprmat)]
dt.meta$Group <- gsub("-", "_", dt.meta$Group)



###############################################
# biological replicates reads pooled together #
#                                             #
#  get exprmat.pool                           #
#                                             #
###############################################

dt.exprSet.pool <- split(dt.meta, by = 'Group')%>%lapply(., function(X){
                         lst.GroupX <- X$ID.sample
                         exprmat.GroupX <- dt.exprmat[,lst.GroupX]
                         exprlist.GroupX <- apply(exprmat.GroupX, 1, mean)
                                
                         Group.ID <- unique(X$Group)
                         dt.GroupX.expr <- data.table(miRNA.ID = names(exprlist.GroupX),
                                                      Group.ID = Group.ID,
                                                      log2CPM = exprlist.GroupX)
                               
                         return(dt.GroupX.expr) 
                                
                         })%>%rbindlist()


dt.exprMat.pool <- dcast(dt.exprSet.pool, miRNA.ID~Group.ID, value.var = 'log2CPM')
exprmat.pool <- data.frame(dt.exprMat.pool, row.names = 1)%>%as.matrix()



###############################################
# get metadata for pooled data                #
#                                             #
#                                             #
###############################################
strsplit(colnames(exprmat.pool), "_")

strsplit(colnames(exprmat.pool), "_")%>%sapply(., '[',1)
dt.meta.pool <- data.table(Group.ID =colnames(exprmat.pool),
                           Organ = strsplit(colnames(exprmat.pool), "_")%>%sapply(., '[',1),
                           Gender = strsplit(colnames(exprmat.pool), "_")%>%sapply(., '[',2),
                           Age = strsplit(colnames(exprmat.pool), "_")%>%sapply(., '[',3))

dt.total.organ.exprSet <- merge(dt.exprSet.pool, dt.meta.pool, by = 'Group.ID')
dt.9.organ.exprSet <- dt.total.organ.exprSet[!Organ%in%c('Tst', 'Utr')]
dt.forPlot <- dcast(dt.9.organ.exprSet, miRNA.ID+Organ+Age~Gender, value.var = 'log2CPM', fun.aggregate = NULL)




###############################################
# colors annotations                          #
#                                             #
# get DEG annotation                                    #
###############################################

dt.gender.DEGs <- readRDS('./RDS/DEG_results/dt.DEG.SexPerOrganAge.cpm.rds')
dt.gender.DEGs$DEGorNot <- ifelse(dt.gender.DEGs$adj.P.Val<0.05&(dt.gender.DEGs$logFC>1 |dt.gender.DEGs$logFC< -1), 'DEG', 'nonDEG')
dt.gender.DEGresult <- data.table(miRNA.ID = dt.gender.DEGs$ID.gene, 
                                  Organ = dt.gender.DEGs$Organ,
                                  Age = dt.gender.DEGs$Age,
                                  DEGtype = dt.gender.DEGs$DEGorNot)

# generate DEG result
dt.gender.DEG.plot <- merge(dt.forPlot, dt.gender.DEGresult, by = c('miRNA.ID', 'Organ', 'Age'))

# Generate DEG list for SupTables
dt.gender.DEG.report <- dt.gender.DEG.plot[DEGtype=='DEG']
write.csv(dt.gender.DEG.report, './supplementary/SupTable_Gender_DEG.csv', row.names = F, quote = F)


# merge DEG result with organ 
dt.gender.DEG.plot$Organ_DEG <- dt.gender.DEG.plot$Organ
dt.gender.DEG.plot[DEGtype=='nonDEG']$Organ_DEG <- 'nonDEG'
dt.gender.DEG.plot$Organ_DEG <- factor(dt.gender.DEG.plot$Organ_DEG, 
                                       levels = c("nonDEG", "Adr", "Brn", "Hrt", "Kdn", "Lng","Lvr", "Msc", "Spl","Thm"))
dt.gender.DEG.plot$DEGtype <- factor(dt.gender.DEG.plot$DEGtype, levels = c('nonDEG', 'DEG'))
dt.gender.DEG.plot$Age <- factor(dt.gender.DEG.plot$Age, levels = c('2', '6', '21', '104'))
dt.gender.DEG.plot$DEG_size <- 1
dt.gender.DEG.plot[DEGtype!='nonDEG']$DEG_size <- 1.5
dt.gender.DEG.plot$DEG_size <- as.numeric(dt.gender.DEG.plot$DEG_size)


# label miRNA.ID
dt.gender.DEG.plot$miRNA.ID.label <- substr(dt.gender.DEG.plot$miRNA.ID, 5, nchar(dt.gender.DEG.plot$miRNA.ID))



# color setting 
color.organ.DEG <- colors.organ[1:9]
color.organ.DEG <- c('gray80', color.organ.DEG)
names(color.organ.DEG) <- c('nonDEG', names(color.organ.DEG)[2:10])







###############################################
# plotting                                    #
#                                             #
#                                             #
###############################################

# four development STAGE


ggplot(dt.gender.DEG.plot, aes(x= F, y= M))+
  geom_point(aes(color=Organ_DEG, size = DEG_size))+
  facet_wrap(~Age, nrow = 2)+
  scale_color_manual(values = color.organ.DEG)+
  theme_bw()+
  scale_x_continuous(limits = c(0, 20))+
  scale_y_continuous(limits = c(0, 20))+
  labs(x= 'Female', y= 'Male')+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=24,face='bold', hjust = 0.5),
        axis.text.y = element_text(size = 24, face = 'bold'),
        axis.title = element_text(size=24,face='bold'),
        legend.title = element_text(size=24,face='bold'),
        legend.text = element_text(size=24,face='bold'),
        strip.text = element_text(size=24,face='bold'))
  
ggsave('./charts/Fig.5b.Scatterplot.gender_r2c2.pdf', width = 16, height = 16)


# nine organs
ggplot(dt.gender.DEG.plot, aes(x= F, y= M))+
  geom_point(aes(color=Organ_DEG, size = DEG_size, shape = Age))+
  facet_wrap(~Organ, nrow = 2)+
  scale_color_manual(values = color.organ.DEG)+
  theme_bw()+
  scale_x_continuous(limits = c(0, 20))+
  scale_y_continuous(limits = c(0, 20))+
  labs(x= 'Female', y= 'Male')+
  theme(axis.text.x = element_text(size=36,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 36, face = 'bold',    color='gray40'),
        axis.title = element_text(size=36,face='bold'),
        legend.title = element_text(size=36,face='bold'),
        legend.text = element_text(size=36,face='bold',color='gray40'),
        strip.text = element_text(size=36,face='bold'))
ggsave('./charts/Fig.5a.Scatterplot.gender.r2c5.pdf', width = 50, height = 20, limitsize = FALSE)




saveRDS(dt.gender.DEG.plot, './RDS/dt.gender.DEG.forPlot.rds')








# test with labels --------------------------------------------------------

dt.gender.DEG.plot$label <- ifelse(dt.gender.DEG.plot$DEGtype=='DEG',
                                   ifelse(dt.gender.DEG.plot$F>10&dt.gender.DEG.plot$M>10 , dt.gender.DEG.plot$miRNA.ID.label, NA), NA)

# dt.gender.DEG.plot$FM_ratio <- dt.gender.DEG.plot$F/ dt.gender.DEG.plot$M
# dt.gender.DEG.plot$label <- ifelse(dt.gender.DEG.plot$FM_ratio>2 | dt.gender.DEG.plot$FM_ratio<0.5, dt.gender.DEG.plot$miRNA.ID, NA )

ggplot(dt.gender.DEG.plot, aes(x= F, y= M))+
  geom_point(aes(color=Organ_DEG, size = DEG_size))+
  facet_wrap(~Age, nrow = 1)+
  scale_color_manual(values = color.organ.DEG)+
  theme_bw()+
  scale_x_continuous(limits = c(0, 20))+
  scale_y_continuous(limits = c(0, 20))+
  labs(x= 'Female', y= 'Male')+
  theme(axis.text.x = element_text(size=36,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 36, face = 'bold',    color='gray40'),
        axis.title = element_text(size=36,face='bold'),
        legend.title = element_text(size=36,face='bold'),
        legend.text = element_text(size=36,face='bold',color='gray40'),
        strip.text = element_text(size=36,face='bold'))
  #geom_text_repel(aes(label = label),box.padding = 1,point.padding = unit(0.4, "lines"), size =5)

ggsave('./charts/Fig.5b.Scatterplot.gender.label.png', width = 42, height = 10, limitsize = FALSE)



ggplot(dt.gender.DEG.plot, aes(x= F, y= M))+
  geom_point(aes(color=Organ_DEG, size = DEG_size, shape = Age))+
  facet_wrap(~Organ, nrow = 2)+
  scale_color_manual(values = color.organ.DEG)+
  theme_bw()+
  scale_x_continuous(limits = c(0, 20))+
  scale_y_continuous(limits = c(0, 20))+
  labs(x= 'Female', y= 'Male')+
  theme(axis.text.x = element_text(size=16,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 16, face = 'bold',    color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=16,face='bold'),
        legend.text = element_text(size=16,face='bold',color='gray40'),
        strip.text = element_text(size=26,face='bold'))+
  geom_text_repel(aes(label = label),box.padding = 1,point.padding = unit(0.4, "lines"), size =5)

ggsave('./charts/Fig.5a.Scatterplot.gender.label.png', width = 30, height = 12, limitsize = FALSE)



# check the expression of miR-200b and miR-429 ----------------------------------------------------

dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
dt.meta <- readRDS('./RDS/metadata.rds')[ID.sample%in%colnames(dt.exprmat)]
dt.meta$Group <- gsub("-", "_", dt.meta$Group)

lst.BrnW2.male.samples <- dt.meta[Organ=='Brn'][Age==2][Sex=='M']$ID.sample
lst.BrnW2.female.samples <- dt.meta[Organ=='Brn'][Age==2][Sex=='F']$ID.sample
lst.TstW2.male.samples <- dt.meta[Organ=='Tst'][Age==2][Sex=='M']$ID.sample


ID <- 'rno-miR-429'
dt.expr.miR.429 <- dt.gender.DEG.plot[miRNA.ID==ID][Organ=='Brn']
dt.expr.miR.429.melt <- melt(dt.expr.miR.429, id.vars = c('Organ', 'Age'), measure.vars = c('F', 'M'), variable.name = 'Sex', value.name = 'logCPM' )

ggplot(dt.expr.miR.429.melt, aes(x= Age, y = logCPM, group = Sex))+
  geom_bar(stat = 'identity', position = 'dodge', aes(fill =Sex ))+
  theme_bw()+
  labs(title = sprintf('Expression of %s', miRNA.ID))+
  theme(plot.title = element_text(size=20,face='bold',color='gray40', hjust = 0.5),
        axis.text.x = element_text(size=18,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 18, face = 'bold',    color='gray40'),
        axis.title = element_text(size=18,face='bold'),
        legend.title = element_text(size=18,face='bold'),
        legend.text = element_text(size=18,face='bold',color='gray40'))
ggsave('./charts/Barplot_miR_429_logCPM.pdf', width = 8, height = 5)  
  
ID <- 'rno-miR-200b-3p'
dt.expr.miR.200 <- dt.gender.DEG.plot[miRNA.ID== ID][Organ=='Brn']
dt.expr.miR.200.melt <- melt(dt.expr.miR.200, id.vars = c('Organ', 'Age'), measure.vars = c('F', 'M'), variable.name = 'Sex', value.name = 'logCPM' )

ggplot(dt.expr.miR.200.melt, aes(x= Age, y = logCPM, group = Sex))+
  geom_bar(stat = 'identity', position = 'dodge', aes(fill =Sex ))+
  theme_bw()+
    labs(title = sprintf('Expression of %s', ID))+
  theme(plot.title = element_text(size=20,face='bold',color='gray40', hjust = 0.5),
        axis.text.x = element_text(size=18,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 18, face = 'bold',    color='gray40'),
        axis.title = element_text(size=18,face='bold'),
        legend.title = element_text(size=18,face='bold'),
        legend.text = element_text(size=18,face='bold',color='gray40'))
ggsave('./charts/Barplot_miR_200b_logCPM.pdf', width = 8, height = 5)  


dt.exprmat.gene <- readRDS('./RDS/exprMat/exprMat_mRNA_log2FPKM.rds')











# novel miRNA  2021.02.07------------------------------------------------------------

dt.novel.miRNA <- readRDS('./RDS/novel_miRNA_analysis/miRNA_info.rds')
lst.novel.miRNA <- dt.novel.miRNA[type.miRNA=='novel']$miRNA.ID

dt.gender.DEG <- dt.gender.DEG.plot[DEGtype == 'DEG']
dt.gender.DEG$type <- ifelse(dt.gender.DEG$F>dt.gender.DEG$M, 'Female', 'Male')
dt.gender.DEG$miRNA.Type <- ifelse(dt.gender.DEG$miRNA.ID%in%lst.novel.miRNA, 'novel', 'annot')
dt.gender.DEG.novel <- dt.gender.DEG[miRNA.Type=='novel']
lst.novel.sex.miRNA <- dt.gender.DEG.novel$miRNA.ID%>%unique()

dt.gender.DEG.novel <- dt.gender.DEG.novel[, .(miRNA.ID, Organ, Age, F, M, type)]%>%unique()
write.csv(dt.gender.DEG.novel, './supplementary/SupTable_Gender_novel.csv',row.names = F)



ggplot(dt.gender.DEG, aes(x= F, y= M))+
  geom_point(aes(color=Organ_DEG, size = DEG_size, shape = miRNA.Type))+
  facet_wrap(~Age, nrow = 1)+
  scale_color_manual(values = color.organ.DEG)+
  theme_bw()+
  scale_x_continuous(limits = c(0, 20))+
  scale_y_continuous(limits = c(0, 20))+
  labs(x= 'Female', y= 'Male')+
  theme(axis.text.x = element_text(size=36,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 36, face = 'bold',    color='gray40'),
        axis.title = element_text(size=36,face='bold'),
        legend.title = element_text(size=36,face='bold'),
        legend.text = element_text(size=36,face='bold',color='gray40'),
        strip.text = element_text(size=36,face='bold'))

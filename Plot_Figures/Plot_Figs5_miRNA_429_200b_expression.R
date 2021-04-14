rm(list = ls())

dt.gender.DEG.plot <- readRDS('./RDS/dt.gender.DEG.forPlot.rds')

# check the expression of miR-200b and miR-429 ----------------------------------------------------

dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
dt.meta <- readRDS('./RDS/metadata.rds')[ID.sample%in%colnames(dt.exprmat)]
dt.meta$Group <- gsub("-", "_", dt.meta$Group)

lst.BrnW2.male.samples <- dt.meta[Organ=='Brn'][Age==2][Sex=='M']$ID.sample
lst.BrnW2.female.samples <- dt.meta[Organ=='Brn'][Age==2][Sex=='F']$ID.sample
lst.TstW2.male.samples <- dt.meta[Organ=='Tst'][Age==2][Sex=='M']$ID.sample


ID <- 'rno-miR-192-5p'
dt.expr.miR.429 <- dt.gender.DEG.plot[miRNA.ID==ID][Organ=='Brn']
dt.expr.miR.429.melt <- melt(dt.expr.miR.429, id.vars = c('Organ', 'Age'), measure.vars = c('F', 'M'), variable.name = 'Sex', value.name = 'logCPM' )

ggplot(dt.expr.miR.429.melt, aes(x= Age, y = logCPM, group = Sex))+
  geom_bar(stat = 'identity', position = 'dodge', aes(fill =Sex ))+
  theme_bw()+
  labs(title = sprintf('Expression of %s', ID))+
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




# check target gene expression --------------------------------------------

dt.exprmat.gene <- readRDS('./RDS/exprMat/exprMat_mRNA_log2FPKM.rds')
dt.miRTarGene <- readRDS('./RDS/dt_MTI_update.rds')
lst.TarGene.200b <- dt.miRTarGene[ID.miRNA == 'rno-miR-200b-3p']$ID.gene
lst.TarGene.429 <- dt.miRTarGene[ID.miRNA == 'rno-miR-429']$ID.gene
lst.TarGene <- union(lst.TarGene.200b, lst.TarGene.429)

ID <- lst.TarGene[1]

dt.exprSet.gene <- melt(dt.exprmat.gene)%>%setnames(., c('ID.gene', 'ID.Sample', 'logFPKM'))%>%as.data.table()
dt.meta.exprSet.gene <- merge(dt.exprSet.gene, dt.meta, by.x = "ID.Sample", by.y = 'ID.sample', all.x = T)
dt.meta.exprSet.gene[, mean.FPKM := mean(logFPKM), by = .(ID.gene, Organ, Sex, Age)]

dt.exprSet <- data.table(ID.gene = dt.meta.exprSet.gene$ID.gene,
                         Organ = dt.meta.exprSet.gene$Organ,
                         Age = dt.meta.exprSet.gene$Age,
                         Sex = dt.meta.exprSet.gene$Sex,
                         logFPKM = dt.meta.exprSet.gene$mean.FPKM)%>%unique()

ID <- 'ENSRNOG00000017863'

dt.forPlot <- dt.exprSet[ID.gene == ID][Organ=='Brn']

ggplot(dt.forPlot, aes(x= Age, y = logFPKM, group = Sex))+
  geom_bar(stat = 'identity', position = 'dodge', aes(fill =Sex ))+
  theme_bw()+
  labs(title = sprintf('Expression of %s', ID))+
  theme(plot.title = element_text(size=20,face='bold',color='gray40', hjust = 0.5),
        axis.text.x = element_text(size=18,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 18, face = 'bold',    color='gray40'),
        axis.title = element_text(size=18,face='bold'),
        legend.title = element_text(size=18,face='bold'),
        legend.text = element_text(size=18,face='bold',color='gray40'))


dt.Gene.name <- readRDS('./RDS/Gene_ENSEMBL_SYMBOL.rds')
dt.Gene.name[SYMBOL=='Zeb1']

dt.miRTarGene[ID.gene== 'ENSRNOG00000017863']

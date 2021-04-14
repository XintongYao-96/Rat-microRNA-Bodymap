library(data.table)
library(magrittr)

dt.cor.allMiRs <- readRDS('./version1_scripts_results/RDS_V1/Corr_5_types/Corr_miRNA_all_miRNA.rds')
dt.cor.sisterMiRs <- readRDS('./version1_scripts_results/RDS_V1/Corr_5_types/Corr_miRNA_sister_miRNA.rds')
dt.cor.host.oppo <- readRDS('./version1_scripts_results/RDS_V1/Corr_5_types/Corr_miRNA_host_gene_oppo.rds')
dt.cor.host.same <- readRDS('./version1_scripts_results/RDS_V1/Corr_5_types/Corr_miRNA_host_gene_same.rds')
dt.cor.targetPair <- readRDS('./version1_scripts_results/RDS_V1/Corr_5_types/Corr_miRNA_target_gene.rds')
dt.cor.allMiRGene <- readRDS( './version1_scripts_results/RDS_V1/Corr_5_types/miRNA_all_gene.rds')


# miRNA and gene correlation based on origin log2CPM (before per organ Z --------

dt.cor.nonZscore <- readRDS('./RDS/correlation_result/correlation_nonZscore.rds')

dt.cor.melt <- melt(dt.cor.nonZscore)%>%setnames(., c('ID.Gene.A', 'ID.Gene.B', 'Corr'))%>%as.data.table()
dt.cor.melt[ , pair := paste0(dt.cor.melt$ID.Gene.A,'_',dt.cor.melt$ID.Gene.B)]


# generate all miRNA and Gene pairs correlation
dt.cor.allmiRGene <- dt.cor.melt
dt.cor.allmiRGene[, Type := 'miRNA-all gene']


# generate miRNA and target gene expression
dt.miRTarGene <- readRDS('./RDS/dt_MTI_update.rds')
pairs.miR.TarGene <- paste0(dt.miRTarGene$ID.miRNA, '_', dt.miRTarGene$ID.gene)

dt.cor.miRTarGene <- dt.cor.melt[pair%in%pairs.miR.TarGene]
dt.cor.miRTarGene[, Type :='miRNA-target gene']


# generate miRNA and host gene (same)
pairs.miRNA.host.same <- paste0(dt.cor.host.same$ID.Gene.A, '_',dt.cor.host.same$ID.Gene.B)%>%unique()
dt.cor.miR.host.same <- dt.cor.melt[pair%in%pairs.miRNA.host.same]
dt.cor.miR.host.same[, Type := 'miRNA-host gene (same)']

# generate miRNA and host gene (oppo)
pairs.miRNA.host.oppo <- paste0(dt.cor.host.oppo$ID.Gene.A, '_', dt.cor.host.oppo$ID.Gene.B)%>%unique()
dt.cor.miR.host.oppo <- dt.cor.melt[pair%in%pairs.miRNA.host.oppo]
dt.cor.miR.host.oppo[, Type := 'miRNA-host gene (opp)']

#miRNA to miRNA CORRELATION
dt.cor.miRNA <- readRDS('./RDS/correlation_result/Corr_miRNA_to_miRNA.rds')
dt.cor.all.miRNA <- dt.cor.miRNA
dt.cor.all.miRNA[, Type := 'miRNA-all miRNA']


# miRNA to sister miRNA
pairs.miRNA.sister <- paste0(dt.cor.sisterMiRs$ID.Gene.A, '_',dt.cor.sisterMiRs$ID.Gene.B)
dt.cor.sister.miRNA <- dt.cor.miRNA[pair%in%pairs.miRNA.sister]
dt.cor.sister.miRNA[, Type := 'miRNA-sister miRNA']

# Plot --------------------------------------------------------------------

dt.forPlot <- rbindlist(list(dt.cor.allmiRGene[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.miRTarGene[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.miR.host.same[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.miR.host.oppo[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.all.miRNA[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.sister.miRNA[, .(ID.Gene.A, ID.Gene.B, Corr, Type)]
                             ))

dt.forPlot$Type <- factor(dt.forPlot$Type,
                          levels=c('miRNA-all gene',
                                   'miRNA-host gene (opp)',
                                   'miRNA-host gene (same)',
                                   'miRNA-target gene',
                                   'miRNA-all miRNA',
                                   'miRNA-sister miRNA')
)

saveRDS(dt.forPlot, './RDS/correlation_result/Corr_allType.rds')


# statistics for text

dt.forPlot[Type=='miRNA-all gene']$Corr%>%median()
dt.forPlot[Type=='miRNA-all gene']$Corr%>%sd()
length(unique(dt.forPlot[Type=='miRNA-host gene (same)']$ID.Gene.A))
length(unique(dt.forPlot[Type=='miRNA-host gene (opp)']$ID.Gene.A))

dt.forPlot[Type=='miRNA-host gene (opp)']$Corr%>%median()
dt.forPlot[Type=='miRNA-host gene (opp)']$Corr%>%sd()

dt.forPlot[Type=='miRNA-host gene (same)']$Corr%>%median()
dt.forPlot[Type=='miRNA-host gene (same)']$Corr%>%sd()

dt.forPlot[Type=='miRNA-sister miRNA']$Corr%>%median()
dt.forPlot[Type=='miRNA-sister miRNA']$Corr%>%sd()



# old version
# p <- ggplot(dt.forPlot,aes(x=Type,y=Corr))+
#   geom_hline(yintercept = 0, linetype='dotted', color='gray60', size=1)+
#   geom_boxplot(outlier.alpha=0,
#                aes(fill=grepl('gene|target',Type)))+
#   labs(x='',y='correlation')+
#   scale_y_continuous(limits=c(-1,1))+
#   theme_bw()+
#   theme(axis.text.x = element_text(size=12,face='bold'),
#         axis.text.y = element_text(size=10, face='bold',color='gray20'),
#         axis.title.x = element_text(size=12,face='bold'))+
#   coord_flip()+
#   guides(fill=FALSE)+
#   scale_fill_manual(values=brewer.pal(3,'Purples')[c(2,1)])
# ggsave('./charts/Fig.old1D.Organ.nonZcore.pdf', p, width=8,height=4)


# 20210216 plot

dt.forPlot <- readRDS('./RDS/correlation_result/Corr_allType.rds')
p <- ggplot(dt.forPlot, aes(x=Type, y = Corr))+
  geom_violin(trim = F, aes(fill=grepl('gene|target',Type), color = grepl('gene|target',Type)))+
  geom_boxplot(outlier.alpha=0, width = 0.1)+
  geom_hline(yintercept = 0, linetype='dotted', color='gray60', size=1)+
  labs(x='',y='correlation')+
  scale_y_continuous(limits=c(-1,1))+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=20,face='bold', hjust = 1, angle = 30),
        axis.text.y = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size=20,face='bold'),
        legend.title = element_text(size=20,face='bold'),
        legend.text = element_text(size=20,face='bold'))+
  guides(fill=FALSE)+
  scale_fill_manual(values=c('#D6A09E','#ADD4EF'))+
  scale_color_manual(values=c('#D6A09E','#ADD4EF'))

ggsave('./charts/Fig.1c.Violin.miRNA_mRNA_cirrelation.pdf', p, width=12,height=6)

#######################################
# plot on per tissue Zscore           #
# miRNA and mRNA expression mat       #
                                      #
#######################################

dt.cor.Zscore <- readRDS('./RDS/correlation_result/correlation_tissue_Z.rds')

# replace the . in rownames with -
lst.miRNA.id <- rownames(dt.cor.Zscore)
lst.miRNA.ID <- gsub('\\.', '-', lst.miRNA.id)
rownames(dt.cor.Zscore) <- lst.miRNA.ID

# melt the correlation matrix
dt.cor.melt <- melt(dt.cor.Zscore)%>%setnames(., c('ID.Gene.A', 'ID.Gene.B', 'Corr'))%>%as.data.table()
dt.cor.melt[ , pair := paste0(dt.cor.melt$ID.Gene.A,'_',dt.cor.melt$ID.Gene.B)]


# generate all miRNA and Gene pairs correlation
dt.cor.allmiRGene <- dt.cor.melt
dt.cor.allmiRGene[, Type := 'miRNA-all gene']


# generate miRNA and target gene expression
dt.miRTarGene <- readRDS('./RDS/dt_MTI_update.rds')
pairs.miR.TarGene <- paste0(dt.miRTarGene$ID.miRNA, '_', dt.miRTarGene$ID.gene)

dt.cor.miRTarGene <- dt.cor.melt[pair%in%pairs.miR.TarGene]
dt.cor.miRTarGene[, Type :='miRNA-target gene']


# generate miRNA and host gene (same)
pairs.miRNA.host.same <- paste0(dt.cor.host.same$ID.Gene.A, '_',dt.cor.host.same$ID.Gene.B)%>%unique()
dt.cor.miR.host.same <- dt.cor.melt[pair%in%pairs.miRNA.host.same]
dt.cor.miR.host.same[, Type := 'miRNA-host gene (same)']

# generate miRNA and host gene (oppo)
pairs.miRNA.host.oppo <- paste0(dt.cor.host.oppo$ID.Gene.A, '_', dt.cor.host.oppo$ID.Gene.B)%>%unique()
dt.cor.miR.host.oppo <- dt.cor.melt[pair%in%pairs.miRNA.host.oppo]
dt.cor.miR.host.oppo[, Type := 'miRNA-host gene (opp)']


dt.cor.miRNA <- readRDS('./RDS/correlation_result/Corr_miRNA_to_miRNA.Z.rds')

# total miRNA
dt.cor.all.miRNA <- dt.cor.miRNA
dt.cor.all.miRNA[, Type := 'miRNA-all miRNA']

# miRNA to sister miRNA
pairs.miRNA.sister <- paste0(dt.cor.sisterMiRs$ID.Gene.A, '_',dt.cor.sisterMiRs$ID.Gene.B)
dt.cor.sister.miRNA <- dt.cor.miRNA[pair%in%pairs.miRNA.sister]
dt.cor.sister.miRNA[, Type := 'miRNA-sister miRNA']

# Plot --------------------------------------------------------------------

dt.forPlot <- rbindlist(list(dt.cor.allmiRGene[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.miRTarGene[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.miR.host.same[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.miR.host.oppo[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.all.miRNA[, .(ID.Gene.A, ID.Gene.B, Corr, Type)],
                             dt.cor.sister.miRNA[, .(ID.Gene.A, ID.Gene.B, Corr, Type)]
))

dt.forPlot$Type <- factor(dt.forPlot$Type,
                          levels=c('miRNA-all gene',
                                   'miRNA-host gene (opp)',
                                   'miRNA-host gene (same)',
                                   'miRNA-target gene',
                                   'miRNA-all miRNA',
                                   'miRNA-sister miRNA')
)


p <- ggplot(dt.forPlot,aes(x=Type,y=Corr))+
  geom_hline(yintercept = 0, linetype='dotted', color='gray60', size=1)+
  geom_boxplot(outlier.alpha=0,
               aes(fill=grepl('gene|target',Type)))+
  labs(x='',y='correlation')+
  scale_y_continuous(limits=c(-1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,face='bold'),
        axis.text.y = element_text(size=10, face='bold',color='gray20'),
        axis.title.x = element_text(size=12,face='bold'))+
  coord_flip()+
  guides(fill=FALSE)+
  scale_fill_manual(values=brewer.pal(3,'Purples')[c(2,1)])

ggsave('./charts/Fig.old1D.OrganZcore.pdf', p, width=8,height=4)

library(data.table) 
library(magrittr)
library(ggplot2)

##############################################################
# chromosome annotation of 725 miRBase annotates miRNAs      #
# 711 from miRBaes gff, 14 searched in NCBI                  #
#                                                            #
# Xintong Yao, 2021.2.10                                     #
##############################################################


# 
# # get miRNA info with chromosome location from miRBase ---------------------------------
# 
# # miRBase_rno.txt is the gff downloaded from miRBase
# dt.miRNA.miRBase.gff <- read.table('./tables/miRBase_rno.txt')
# lst.ID <- dt.miRNA.miRBase.gff$V9
# miRNA.ID<- strsplit(lst.ID, ';')%>%sapply(., '[',3)%>%substr(., 6, nchar(.))
# dt.miRNA.miRBase.gff$ID <- miRNA.ID
# 
# dt.miRNA.miRBase <- data.table(miRNA.name = dt.miRNA.miRBase.gff$ID,
#                                   miRNA.type = dt.miRNA.miRBase.gff$V3,
#                                   chromosome.name = dt.miRNA.miRBase.gff$V1)%>%unique() # 1275 miRNA Annotation including pri-miRNA
# 
# miRNA.number <- table(dt.miRNA.miRBase$miRNA.name)
# miRNA.number[miRNA.number>1]
# 
# # get miRNAs includede in our exprmat
# dt.miRNA.info <- readRDS('./RDS/novel_miRNA_analysis/miRNA_info.rds')
# dt.miRNA.exprmat <- data.table(miRNA.name = dt.miRNA.info[type.miRNA=='annotated']$miRNA.ID)
# 
# dt.miRNA.annot <- merge(dt.miRNA.exprmat, dt.miRNA.miRBase, by = 'miRNA.name', all.x = T)
# # miRNA.number <- table(dt.miRNA.annot$miRNA.name)
# # miRNA.number[miRNA.number>1]
# 
# write.csv(dt.miRNA.annot, './tables/miRNA_chromosome_annotaion_forEdit.csv',row.names = F)
# 
# 
# 
# # annotation of total miRNA -----------------------------------------------
# 
# dt.miRNA.annotation <-fread('./tables/miRNA_chromosome_annotaion_Edit.csv')
# dt.miRNA.annotation$type.miRNA <- 'annotated'
# 
# dt.miRNA.info$chromosome.name <- strsplit(dt.miRNA.info$miRNA.ID, '-')%>%sapply(., '[',2)
# colnames(dt.miRNA.info)[1] <- 'miRNA.name'
# dt.miRNA.total.annotation <- rbind(dt.miRNA.info[type.miRNA=='novel'], 
#                                    dt.miRNA.annotation[, .(miRNA.name, type.miRNA, chromosome.name)])
# saveRDS(dt.miRNA.total.annotation, './RDS/miRNA_chromosome_annotation_total.rds')
# 



# Plot --------------------------------------------------------------------

dt.miRNA.total.annotation <-  readRDS('./RDS/miRNA_Chromosome_annot.rds') 

dt.N.miRNA <- dt.miRNA.total.annotation[, .(N.miRNA = .N), by = .(type.miRNA, chromosome.name)]


dt.N.miRNA <- dt.N.miRNA[chromosome.name!= 'NA']
dt.N.miRNA$chromosome.name <- factor(dt.N.miRNA$chromosome.name, 
                                     levels = paste0('chr', c(seq(1,20), 'X')))
dt.N.miRNA$type.miRNA <- factor(dt.N.miRNA$type.miRNA, 
                                levels = rev(c('annotated','novel')))

ggplot(dt.N.miRNA, aes(x = chromosome.name, y = N.miRNA, group = type.miRNA))+
  geom_bar(stat = 'identity', position="stack", aes(fill=type.miRNA ))+
  labs(x= 'Chromosome', y = 'Number of miRNAs detected')+
  theme_bw()+
  scale_fill_manual(values = c('#a2d5f2','#07689f'))+
  theme(axis.text.x = element_text(size=16,face='bold',color='black', angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 16, face = 'bold'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=16,face='bold'),
        legend.text = element_text(size=16,face='bold'))
ggsave('./charts/Fig7b_Barplot_number_miRNA_chromosome.pdf', width = 14, height = 5)

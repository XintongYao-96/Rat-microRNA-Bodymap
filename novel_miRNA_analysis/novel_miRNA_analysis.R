library(data.table)
library(magrittr)


dt.exprmat <- readRDS('./RDS/exprMat/exprMat_miRNA_log2CPM.rds')
lst.miRNA <- rownames(dt.exprmat)


# generate novel miRNA
lst.novel.miRNA <- lst.miRNA[grep('rno-chr', lst.miRNA)]
lst.chr <- strsplit(lst.novel.miRNA, '-')%>%sapply(., '[',2)

lst.miRNA.type <- ifelse(lst.miRNA%in%lst.novel.miRNA, 'novel', 'annotated')
dt.miRNA.info <- data.table(miRNA.ID = lst.miRNA,
                            type.miRNA = lst.miRNA.type)
saveRDS(dt.miRNA.info, './RDS/novel_miRNA_analysis/miRNA_info.rds')


# PLOT NOVEL miRNA Number on each chromosome
dt.novel.miRNA.info <- data.table(miRNA.ID = lst.novel.miRNA,
                                  chr = lst.chr)
dt.novel.miRNA.info$chr.ID <- substr(as.character(dt.novel.miRNA.info$chr), 4, 5)
dt.N.miRNA.chromo <- dt.novel.miRNA.info[, .(N.miRNA = .N), by = 'chr.ID']
dt.N.miRNA.chromo$chr.ID <- factor(dt.N.miRNA.chromo$chr.ID ,
                                   levels = c(seq(1:20), 'X'))

ggplot(dt.N.miRNA.chromo, aes(x = chr.ID, y = N.miRNA))+
  geom_bar(stat = 'identity', fill = '#5C6994')+
  labs(x= 'Chromosome', y = 'Number of novel miRNA')+
  theme_bw()+
  theme(axis.text.x = element_text(size=16,face='bold',color='gray40', hjust = 0.5),
        axis.text.y = element_text(size = 16, face = 'bold'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=16,face='bold'),
        legend.text = element_text(size=16,face='bold'))
ggsave('./charts/Novel_miRNA_number_chromosome.pdf', width = 8, height = 4)  



library(data.table)
library(magrittr)
library(reshape2)

dt.miRNA.z.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_miRNA_log2CPM.rds')
dt.mRNA.z.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_mRNA_log2FPKM.rds')
dt.cor <- readRDS('./RDS/correlation_result/cor_all_miRNA_mRNA.rds')%>%t()
dt.meta <- readRDS('./RDS/metadata.rds')[ID.sample%in%colnames(dt.miRNA.z.exprmat)]

dt.cor.melt <- melt(dt.cor)%>%setnames(., c('gene', 'miRNA', 'cor'))%>%as.data.table()

name.gene <- 'ENSRNOG00000031090'
name.miRNA <- 'rno.miR.29a.3p'

dt.cor[name.gene,name.miRNA]
dt.cor.melt[gene==name.gene][miRNA==name.miRNA]


dt.order.cor <- dt.cor.melt[order(cor)]



dt.neg.cor <- dt.order.cor[1:20]

dt.neg.cor$gene <- dt.neg.cor$gene%>%as.character(.)
dt.neg.cor$miRNA <- dt.neg.cor$miRNA%>%as.character(.)


# corrlation
split(dt.neg.cor, by=c('gene', 'miRNA'))%>%lapply(., function(x){
  name.gene <- x$gene
  name.gene
  name.miRNA <- x$miRNA
  name.miRNA
  
  lst.gene <- dt.mRNA.z.exprmat[name.gene,]
  lst.miRNA <- dt.miRNA.z.exprmat[name.miRNA,]
  
  
  dt.forPlot <- data.table(miRNA=lst.miRNA,
                           gene=lst.gene,
                           ID.sample = names(lst.gene))
  dt.forPlot <- merge(dt.forPlot, dt.meta, by = 'ID.sample')
  
  ggplot(dt.forPlot, aes(x=miRNA, y=gene))+
    geom_point(aes(color = Age))+
    theme_bw()+
    labs(title = sprintf('cor = %.3f%%', x$cor))+
    theme(axis.text = element_text(size=16,face='bold',color='gray40'),
          axis.title = element_text(size=16,face='bold'),
          legend.title = element_text(size=16,face='bold'),
          legend.text = element_text(size=16,face='bold',color='gray40'))
  
  dir = paste0(c('./charts/correlation_scatterplots/'),'neg',name.gene,"_",name.miRNA,".png")
  ggsave(dir, width = 6, height = 5)
  
  
})



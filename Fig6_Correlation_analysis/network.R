library(data.table)
library(magrittr)

dt.miRTarGene <- readRDS('./RDS/dt_MTI_update.rds')
lst.cluster1.miRNA <- readRDS('./RDS/correlation_result/tier1.miRNA.rds')
lst.cluster2.miRNA <- readRDS('./RDS/correlation_result/tier2.miRNA.rds')
lst.cluster1.gene <- readRDS('./RDS/correlation_result/tier1.mRNA.rds')
lst.cluster2.gene <- readRDS('./RDS/correlation_result/tier2.mRNA.rds')



dt.cor <- readRDS('./RDS/correlation_result/correlation_tissue_Z.rds')
dt.melt.cor <- melt(dt.cor)%>%setnames(., c('ID.miRNA' ,'ID.gene','cor'))
dt.melt.cor$ID.miRNA <- gsub('\\.', '-', dt.melt.cor$ID.miRNA)

#merge(dt.cluster1.miRTarGene, dt.melt.cor)


##########################################
#  generate miRNA-gene for visualization #
#                                        #
##########################################

lst.cluster.miRNA  <- union(lst.cluster1.miRNA, lst.cluster2.miRNA)
delete.miRNA <- which(lst.cluster.miRNA=="rno-miR-125b-5p")
lst.cluster.miRNA <- lst.cluster.miRNA[-delete.miRNA]

lst.cluster.gene <- union(lst.cluster1.gene, lst.cluster2.gene)


dt.cluster.miRTarGene <- dt.miRTarGene[ID.miRNA%in%lst.cluster.miRNA][ID.gene%in%lst.cluster.gene]
dt.cluster.miRTarGene.cor <- merge(dt.cluster.miRTarGene, dt.melt.cor)
dt.cluster.miRTarGene.cor <- dt.cluster.miRTarGene.cor[,.(ID.miRNA, ID.gene, cor)]

write.csv(dt.cluster.miRTarGene.cor, './tables/dt.cluster.miRTarGene.cor.csv', row.names = F, quote = F)


dt.cluster.miRTarGene.Neg.cor <- dt.cluster.miRTarGene.cor[cor<0]





# pathway -----------------------------------------------------------------

Pathway.tier1.gene <- readRDS('./RDS/Pathway_result/Pathway.tier1.gene.rds')
Pathway.tier2.gene <- readRDS('./RDS/Pathway_result/Pathway.tier2.gene.rds')

View(pathway.tier1.gene)

SupTable.pathway.tier1.gene <- pathway.tier1.gene[1:30, c('ID', 'Description', 'p.adjust')]
write.csv(SupTable.pathway.tier1.gene,'./supplementary/SupTable.pathway.cluster1.gene.csv', row.names = F, quote = F)

SupTable.pathway.tier2.gene <- pathway.tier2.gene[1:30, c('ID', 'Description', 'p.adjust')]
write.csv(SupTable.pathway.tier2.gene,'./supplementary/SupTable.pathway.cluster2.gene.csv', row.names = F, quote = F)



######################################
#  select pathway for visualization  #
#                                    #
######################################

# tier 1 pathway

order <- seq(1,10)
dt.gene.in.tier1pathway <- sapply(order, function(x){
                              name.pathway <- Pathway.tier1.gene[x,]$Description%>%unique(.)
                              lst.gene.pathway <- unlist(strsplit(Pathway.tier1.gene$geneID[x], split="/"))
          
                              lst.TaGene.pathway <- intersect(dt.cluster.miRTarGene.Neg.cor$ID.gene,lst.gene.pathway )
                              
                              dt.gene.in.pathway = data.table(pathway = name.pathway,
                                                              ID.gene = lst.TaGene.pathway)
                              
                              return(list(dt.gene.in.pathway))
                              
                            })%>%rbindlist()

length(unique(dt.gene.in.tier1pathway$ID.gene))

order <- seq(1,18)
dt.gene.in.tier2pathway <- sapply(order, function(x){
  name.pathway <- Pathway.tier2.gene[x,]$Description%>%unique(.)
  lst.gene.pathway <- unlist(strsplit(Pathway.tier2.gene$geneID[x], split="/"))
  
  lst.TaGene.pathway <- intersect(dt.cluster.miRTarGene.Neg.cor$ID.gene,lst.gene.pathway )
  
  dt.gene.in.pathway = data.table(pathway = name.pathway,
                                  ID.gene = lst.TaGene.pathway)
  
  return(list(dt.gene.in.pathway))
  
})%>%rbindlist()
# 
# lst.gene.antigen <- unlist(strsplit(Pathway.tier2.gene$geneID[1], split="/"))
# intersect(dt.cluster.miRTarGene.cor$ID.gene, lst.gene.antigen)

dt.gene.pathway <- rbind(dt.gene.in.tier1pathway,dt.gene.in.tier2pathway)

saveRDS(dt.gene.in.tier1pathway, './RDS/Pathway_result/dt.pathway.cluster1.gene.rds')
saveRDS(dt.gene.in.tier2pathway, './RDS/Pathway_result/dt.pathway.cluster2.gene.rds')



###########################################
#  combine miRNA-gene-pathway interaction #
#        for Cytosacape                   #
###########################################

# generate network.csv
dt.miRNA.gene <- data.table(source = dt.cluster.miRTarGene.Neg.cor$ID.miRNA,
                            target = dt.cluster.miRTarGene.Neg.cor$ID.gene,
                            value =  dt.cluster.miRTarGene.Neg.cor$cor)

dt.pathway.gene <- data.table(source = dt.gene.pathway$pathway,
                              target =dt.gene.pathway$ID.gene,
                              value = 1)
dt.network <- rbind(dt.miRNA.gene, dt.pathway.gene)

write.csv(dt.network, './tables/dt.network.csv', row.names = F, quote = F)
saveRDS(dt.pathway.gene, './RDS/Pathway_result/dt.pathway.gene.rds')


# generate node.csv
dt.miRNA.node <- data.table(node = dt.miRNA.gene$source,
                            type = 'miRNA')
dt.gene.node <-  data.table(node = dt.miRNA.gene$target,
                            type = 'gene')
dt.pathway.node <-  data.table(node = dt.pathway.gene$source,
                            type = 'pathway')

dt.node.attribute <- rbind(dt.miRNA.node,dt.gene.node, dt.pathway.node )

write.csv(dt.node.attribute, './tables/dt.node.attribute.csv', quote = F)

library(data.table)
library(magrittr)

dt.miRTarGene <- readRDS('./RDS/dt_MTI_update.rds')

lst.cluster1.miRNA <- readRDS('./RDS/correlation_result/tier1.miRNA.rds')
lst.cluster2.miRNA <- readRDS('./RDS/correlation_result/tier2.miRNA.rds')
lst.cluster1.gene <- readRDS('./RDS/correlation_result/tier1.mRNA.rds')
lst.cluster2.gene <- readRDS('./RDS/correlation_result/tier2.mRNA.rds')


lst.cluster1.miRNA <- gsub('\\.', '-', lst.cluster1.miRNA)
lst.cluster2.miRNA <- gsub('\\.', '-', lst.cluster2.miRNA)


# get target gene for cluster1&2miRNA

dt.cluster1.miRTarGene <- dt.miRTarGene[ID.miRNA%in%lst.cluster1.miRNA]
dt.cluster2.miRTarGene <- dt.miRTarGene[ID.miRNA%in%lst.cluster2.miRNA]
lst.cluster1.miRTarGene <- dt.cluster1.miRTarGene$ID.gene%>%unique()
lst.cluster2.miRTarGene <- dt.cluster2.miRTarGene$ID.gene%>%unique()

# intersect of target gene and negative correlation gene
lst.cluster2gene.Tar.NegCor <- intersect(lst.cluster2.miRTarGene, lst.cluster2.gene)
lst.cluster1gene.Tar.NegCor <- intersect(lst.cluster1.miRTarGene, lst.cluster1.gene)



# pathway enrichment analysis ---------------------------------------------

library(clusterProfiler)
library(org.Rn.eg.db)

Pathway.tier1.gene <- data.frame(enrichGO(lst.cluster1.gene, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))
Pathway.tier2.gene <- data.frame(enrichGO(lst.cluster2.gene, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))
Pathway.cluster1gene.Tar.NegCor <- data.frame(enrichGO(lst.cluster1gene.Tar.NegCor, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))
Pathway.cluster2gene.Tar.NegCor <- data.frame(enrichGO(lst.cluster2gene.Tar.NegCor, 'org.Rn.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.05, pAdjustMethod = "fdr", qvalueCutoff = 1))

saveRDS(Pathway.tier1.gene, './RDS/Pathway_result/Pathway.tier1.gene.rds')
saveRDS(Pathway.tier2.gene, './RDS/Pathway_result/Pathway.tier2.gene.rds')


lst.gene.chromosome.segregation <- unlist(strsplit(Pathway.tier1.gene$geneID[1], split="/"))
lst.gene.chromosome.segregation.withTarget <- intersect(lst.gene.chromosome.segregation, lst.cluster1gene.Tar.NegCor)

dt.cluster1.miRTarGene[ID.gene%in%lst.gene.chromosome.segregation.withTarget]

write.csv(Pathway.tier1.gene, './tables/cluster1Gene_Pathway.csv', quote = F, row.names = F)
write.csv(Pathway.tier2.gene, './tables/cluster2Gene_Pathway.csv', quote = F, row.names = F)
write.csv(Pathway.cluster1gene.Tar.NegCor, './tables/cluster1miRNA_targetGene_NegCor_Pathway.csv', quote = F, row.names = F)
write.csv(Pathway.cluster2gene.Tar.NegCor, './tables/cluster2miRNA_targetGene_NegCor_Pathway.csv', quote = F, row.names = F)



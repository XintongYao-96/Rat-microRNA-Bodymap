rm(list = ls())
library(data.table)
library(magrittr)

dt.exprSet <- readRDS('./RDS/exprMat/exprSet_rawcount_1033miR320Sample.rds')
dt.metadata <- readRDS('./RDS/metadata.rds')


# Number of miRNA detected for each sample --------------------------------

dt.N.total.miRNA <- dt.exprSet[Count >3, .(Total.miRNA.Number = .N), by = .(Sample.ID)]

IDs.novel.miRNA <- unique(dt.exprSet$miRNA.ID)%>%grep('chr', ., value = T)
IDs.known.miRNA <- setdiff(unique(dt.exprSet$miRNA.ID), IDs.novel.miRNA)

dt.exprSet.known <- dt.exprSet[miRNA.ID%in%IDs.known.miRNA]
dt.N.known.miRNA <- dt.exprSet.known[Count >3, .(Known.miRNA.Number = .N), by = .(Sample.ID)]

dt.exprSet.novel <- dt.exprSet[miRNA.ID%in%IDs.novel.miRNA]
dt.N.novel.miRNA <- dt.exprSet.novel[Count >3, .(Novel.miRNA.Number = .N), by = .(Sample.ID)]

dt.meta.merge1 <- merge(dt.metadata, dt.N.total.miRNA, by.x = 'Colnames.miRNA', by.y = 'Sample.ID' )
dt.meta.merge2 <- merge(dt.meta.merge1, dt.N.known.miRNA, by.x = 'Colnames.miRNA', by.y = 'Sample.ID' )
dt.meta.merge <- merge(dt.meta.merge2, dt.N.novel.miRNA, by.x = 'Colnames.miRNA', by.y = 'Sample.ID' )

saveRDS(dt.meta.merge, './RDS/metadata_with_miRNA_Number.rds')



# Number of miRNA detected for each chromosome ----------------------------

dt.chromosome.annot <- readRDS('./RDS/miRNA_chromosome_annotation_total.rds')
dt.known.miR.chromosome.annot <- dt.chromosome.annot[type.miRNA=="annotated" ]

chroms.novel.miRNA <- lapply(strsplit(IDs.novel.miRNA, '\\_'),  '[', 3)%>%as.character()

dt.novel.miR.chromosome.annot <- data.table(miRNA.name = IDs.novel.miRNA,
                                            type.miRNA = 'novel',
                                            chromosome.name = chroms.novel.miRNA)

dt.miR.chromosome.annot <- rbind(dt.novel.miR.chromosome.annot, dt.known.miR.chromosome.annot)

saveRDS(dt.miR.chromosome.annot, './RDS/miRNA_Chromosome_annot.rds')


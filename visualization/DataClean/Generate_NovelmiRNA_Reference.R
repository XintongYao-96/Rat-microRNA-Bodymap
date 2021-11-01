#######################################
# Generate reference of novel miRNA   #
# fasta file                          #
#                                     #
# Xintong Yao                         #
#######################################

rm(list = ls())
library(data.table)
library(magrittr)

dt.miRDeep.result <- fread('./tables/miRDeepTable_NovelMiRNA_merge320.csv')

# Filtering:
# Criteria: score >5 and p-value<0.05 
dt.pass.novel.miRNA <- dt.miRDeep.result[`miRDeep2 score`>=4][`significant randfold p-value`=='yes'][`rfam alert`!= 'rRNA']
IDs.novel.miRNA <- paste0('rno_miR_', dt.pass.novel.miRNA$`provisional id`)
IDs.novel.premiRNA <- paste0('rno_mir_', dt.pass.novel.miRNA$`provisional id`)

dt.novel.info <- data.table(miRNA_ID = IDs.novel.miRNA,
                            mature_seq = toupper(dt.pass.novel.miRNA$`consensus mature sequence`),
                            star_seq = toupper(dt.pass.novel.miRNA$`consensus star sequence`),
                            pre_seq = toupper(dt.pass.novel.miRNA$`consensus precursor sequence`),
                            position = toupper(dt.pass.novel.miRNA$`precursor coordinate`)
                            )
saveRDS(dt.novel.info, './RDS/metadata_novel_miRNA.RDS')
write.csv(dt.novel.info, './tables/Table2_Metadata_Novel_miRNA_Information.csv', quote = F, row.names = F)


# Write fasta -------------------------------------------------------------

SEQs.novel.miRNA <- gsub('U', 'T', toupper(dt.pass.novel.miRNA$`consensus mature sequence`))
SEQs.novel.pre.miRNA <- gsub('U', 'T', toupper(dt.pass.novel.miRNA$`consensus precursor sequence`))


N.rows = 2*nrow(dt.novel.info)
Index.odd <- seq(1, N.rows-1, 2)
Index.even <- seq(2, N.rows, 2)

# mature miRNA
Lst.reference <- rep(0, N.rows)
Lst.reference[Index.odd]<- paste0('>', IDs.novel.miRNA)
Lst.reference[Index.even]<- SEQs.novel.miRNA
head(Lst.reference)
write.table(Lst.reference, './tables/rno_mature_novel.fa',sep = '\t', quote = F, row.names = F)

# precursor miRNA
Lst.premiRNA.ref <- rep(0, N.rows)
Lst.premiRNA.ref[Index.odd]<- paste0('>', IDs.novel.premiRNA)
Lst.premiRNA.ref[Index.even]<- SEQs.novel.pre.miRNA
write.table(Lst.premiRNA.ref, './tables/rno_precursor_novel.fa',sep = '\t', quote = F, row.names = F)


# validation --------------------------------------------------------------

grep('chrU', IDs.novel.miRNA)
grep('chrUn|chrM',IDs.novel.miRNA)

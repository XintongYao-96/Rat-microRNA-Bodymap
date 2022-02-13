library(data.table)
library(magrittr)

dt.exprset.novel <- readRDS('./RDS/exprMat/exprSet_novel_count_2022.rds')
dt.exprset.known <- readRDS('./RDS/exprMat/exprSet_2022.rds')

dt.exprmat.novel <- fread('./tables/Table2_Meta_Novel_miRNA_Information_2022.csv')
dt.exprmat.known <- readRDS('./RDS/exprMat/exprMat_CPM_r604c318.rds')
dt.meta <- fread('./tables/Table1_metadata_forPublish_220207.csv')
dt.meta[, Colnames_exprmat_pub := paste0(dt.meta$Sample_ID, '_', dt.meta$SRA_ID)]

IDs.know.include <- rownames(dt.exprmat.known)
IDs.novel.include <- dt.exprmat.novel$miRNA_ID

dt.exprset.combine <- rbind(dt.exprset.novel[miRNA.ID%in%IDs.novel.include],
                            dt.exprset.known[miRNA.ID%in%IDs.know.include])

dt.exprmat.combine <- dcast(dt.exprset.combine, miRNA.ID~Sample.ID, value.var = 'Count')
dt.exprmat <- data.frame(dt.exprmat.combine, row.names = 1)%>%as.matrix()%>%replace(., is.na(.),0)



# Set the order for the expression matrix ---------------------------------

# Order of miRNA: Decreasing mean expression
mean_expr <- apply(dt.exprmat, 1, mean)
order.miRNA <- names(sort(-mean_expr))

# Order of sample: Increasing Organ, Age and sex
dt.meta.order <- dt.meta[order(Organ, Age, Sex)]
order.sample <- dt.meta.order$Colnames.miRNA

dt.exprmat.order <- dt.exprmat[order.miRNA,order.sample]
colnames(dt.exprmat.order) <- dt.meta.order$Colnames_exprmat_pub



# Annot sequence ----------------------------------------------------------

dt.seq.known <- fread('./reference/miRNA_expand/mature_align2_hairpin.txt')

dt.seq.known.forCombine <- data.table(miRNA.ID = dt.seq.known$V1,
                                      seq = gsub('T', 'U', dt.seq.known$V10))
dt.seq.known.forCombine <- dt.seq.known.forCombine[miRNA.ID%in%order.miRNA]
dt.seq.novel.forCombine <- data.table(miRNA.ID = dt.exprmat.novel$miRNA_ID,
                                      seq = dt.exprmat.novel$mature_seq)
dt.seq.combined <- rbind(dt.seq.known.forCombine, dt.seq.novel.forCombine)

df.exprmat.order <- data.frame(dt.exprmat.order)
df.exprmat.order$miRNA.ID <- rownames(df.exprmat.order)
df.exprmat.merge <- merge(df.exprmat.order, dt.seq.combined, by = 'miRNA.ID', all.x = T)

df.exprmat.merge <- df.exprmat.merge[, c(colnames(df.exprmat.merge)[322], colnames(df.exprmat.merge)[1:321])]
dt.exprmat.merge <- data.table(df.exprmat.merge)
setkey(dt.exprmat.merge, 'miRNA.ID')
dt.exprmat <- dt.exprmat.merge[order.miRNA]

saveRDS(dt.exprmat, './RDS/exprmat_forPublish.rds')
write.csv(dt.exprmat, './tables/Exprmat_forPublish.csv', row.names = F, quote = F)



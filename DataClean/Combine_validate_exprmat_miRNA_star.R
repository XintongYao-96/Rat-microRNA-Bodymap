library(data.table)
library(magrittr)

dt.exprmat.valid <- fread('./tables/malerat55_miRNA_filtered_unimean_r424c180.txt')
dt.exprmat.valid <- data.frame(dt.exprmat.valid, row.names = 1)%>%as.matrix()


dt.exprset.valid <- melt(dt.exprmat.valid)%>%setnames(., c('miRNA.ID', 'Sample.ID', 'logI'))%>%as.data.table()
lst.miRNA.ID <- dt.exprset.valid$miRNA.ID
dt.exprset.valid$group.ID <- dt.exprset.valid$miRNA.ID
dt.exprset.valid[grep('\\*$',lst.miRNA.ID)]$group.ID <- gsub('\\*', '', dt.exprset.valid[grep('\\*$',lst.miRNA.ID)]$miRNA.ID)
dt.exprset.valid$logI <- as.numeric(dt.exprset.valid$logI)


dt.exprset.valid.combine <- dt.exprset.valid[, .(Intensity = mean(logI)), by = .(group.ID, Sample.ID)]
dt.exprmat.valid.combine <- dcast(dt.exprset.valid.combine, group.ID~Sample.ID, value.var='Intensity')

exprmat.valid.combine <- data.frame(dt.exprmat.valid.combine, row.names = 1)%>%as.matrix()

saveRDS(exprmat.valid.combine, './RDS/exprMat/exprmat_validate_combine_miRNA_star.rds')

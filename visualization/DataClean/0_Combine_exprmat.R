rm(list = ls())
library(data.table)
library(magrittr)

dir.readcount <- './quantification/'
folders.readcount <- dir(dir.readcount)


dt.exprSet <- lapply(folders.readcount, function(x){
  
  dir.currFile <- sprintf('%s/%s/miRNAs_expressed_all_samples_%s.csv', dir.readcount, x,x)
  dt.in <- fread(dir.currFile)
  dt.in$Sample.ID <-x
  
  return(dt.in)
})%>%rbindlist()


dt.count.exprSet <- data.table(miRNA.ID = dt.exprSet$`#miRNA`,
                               Sample.ID = dt.exprSet$Sample.ID,
                               count = dt.exprSet$read_count)

dt.count.exprSet.unique <- dt.count.exprSet[, .(Count = max(count)), by = .(miRNA.ID, Sample.ID)]
dt.count.exprSet.unique.dcast <- dcast(dt.count.exprSet.unique, miRNA.ID~Sample.ID, fill = 'Count')
mat.count.exprMat <- data.frame(dt.count.exprSet.unique.dcast, row.names = 1)%>%as.matrix()

saveRDS(dt.count.exprSet.unique, './RDS/exprMat/exprSet_rawcount_1033miR320Sample.rds')
saveRDS(mat.count.exprMat, './RDS/exprMat/exprMat_rawcount_r1033c320.rds')
write.csv(mat.count.exprMat, './tables/exprMat_rawcount_r1033c320.csv', row.names = T, quote = F)


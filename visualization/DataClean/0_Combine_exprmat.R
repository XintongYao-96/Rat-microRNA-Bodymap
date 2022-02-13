rm(list = ls())
library(data.table)
library(magrittr)


#  miRNA (annotated) ------------------------------------------------------

dir.readcount <- './data_202201/miRNA_readcount'
folders.readcount <- dir(dir.readcount)
dt.meta <- readRDS('./RDS/metadata.rds')

Name.Sample <- dt.meta$Colnames.miRNA

dt.exprSet <- lapply(Name.Sample, function(x){
  
  dir.currFile <- sprintf('%s/%s.miRBase.ReadCount', dir.readcount, x,x)
  dt.in <- fread(dir.currFile)
  dt.in$Sample.ID <-x
  colnames(dt.in)[1:2] <- c('miRNA.ID', 'Count')
  
  return(dt.in)
})%>%rbindlist()


dt.exprMat <- dcast.data.table(dt.exprSet, miRNA.ID~Sample.ID, value.var  = 'Count')
exprMat <- data.frame(dt.exprMat, row.names = 1)%>%as.matrix()%>%replace(., is.na(.),0)


# Get the colnames free of lane ID
sample.order <- dt.meta$Colnames.miRNA
exprmat.order <- exprMat[, sample.order]
colnames(exprmat.order) <-  dt.meta$ID.sample

saveRDS(dt.exprSet, './RDS/exprMat/exprSet_2022.rds')
saveRDS(exprmat.order, './RDS/exprMat/exprMat_rawcount_r745c320.rds')





# readstats ---------------------------------------------------------------

dir.readstats <- './data_202201/readstats'
folders.readstats <- dir(dir.readstats)

Name.Sample <- dt.meta$Colnames.miRNA
# x <- Name.Sample[1]
dt.ReadStats <- lapply(Name.Sample, function(x){
  
  dir.currFile <- sprintf('%s/%s.ReadStats', dir.readstats, x)
  dt.in <- fread(dir.currFile)
  dt.in$Sample.ID <-x
  colnames(dt.in)[1:2] <- c('Type', 'Count')
  
  return(dt.in)
})%>%rbindlist()


saveRDS(dt.ReadStats, './RDS_202201/ReadStats.rds')


# Novel miRNA  ------------------------------------------------------------

dir.readcount <- './data_202201/novel_readcount'
folders.readcount <- dir(dir.readcount)
dt.meta <- readRDS('./RDS/metadata.rds')

Name.Sample <- dt.meta$Colnames.miRNA

dt.exprSet <- lapply(Name.Sample, function(x){
  
  dir.currFile <- sprintf('%s/%s.novel_miRNA.ReadCount', dir.readcount,x)
  dt.in <- fread(dir.currFile)
  dt.in$Sample.ID <-x
  colnames(dt.in)[1:2] <- c('miRNA.ID', 'Count')
  
  return(dt.in)
})%>%rbindlist()


dt.exprMat <- dcast.data.table(dt.exprSet, miRNA.ID~Sample.ID, value.var  = 'Count')
exprMat <- data.frame(dt.exprMat, row.names = 1)%>%as.matrix()%>%replace(., is.na(.),0)


# Get the colnames free of lane ID
sample.order <- dt.meta$Colnames.miRNA
exprmat.order <- exprMat[, sample.order]
colnames(exprmat.order) <-  dt.meta$ID.sample

saveRDS(exprmat.order, './RDS/exprMat/exprMat_novel_count_r270c320.rds')

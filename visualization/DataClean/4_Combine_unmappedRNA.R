rm(list = ls())
library(data.table)
library(magrittr)


#  miRNA (annotated) ------------------------------------------------------

dir.readcount <- './data_202201/RNAunmapped_readcount'
folders.readcount <- dir(dir.readcount)
dt.meta <- readRDS('./RDS/metadata.rds')

Name.Sample <- dt.meta$Colnames.miRNA

dt.exprSet <- lapply(Name.Sample, function(x){
  dir.currFile <- sprintf('%s/%s.TranscriptUnmapped.ReadCount', dir.readcount, x)
  dt.in <- fread(dir.currFile)[V2>100]
  dt.in$Sample.ID <-x
  colnames(dt.in)[1:2] <- c('Sequence', 'Count')
  
  return(dt.in)
})%>%rbindlist()


dt.exprMat <- dcast.data.table(dt.exprSet, Sequence~Sample.ID, value.var  = 'Count')
exprMat <- data.frame(dt.exprMat, row.names = 1)%>%as.matrix()

# Get the colnames free of lane ID
sample.order <- dt.meta$Colnames.miRNA
exprmat.order <- exprMat[, sample.order]
colnames(exprmat.order) <-  dt.meta$ID.sample

saveRDS(dt.exprSet, './RDS/exprMat/exprSet_RNAunmapped.rds')
saveRDS(exprmat.order, './RDS/exprMat/exprMat_RNAunmapped_rawcount_r3393c320.rds')



# Get fasta files
N.rows = 2*nrow(exprmat.order)
Index.odd <- seq(1, N.rows-1, 2)
Index.even <- seq(2, N.rows, 2)

# mature miRNA
Lst.reference <- rep(0, N.rows)
Lst.reference[Index.odd]<- paste0('>', rownames(exprmat.order))
Lst.reference[Index.even]<- rownames(exprmat.order)
head(Lst.reference)
write.table(Lst.reference, './tables/RNAunmaped.fa',sep = '\t', quote = F, row.names = F, col.names = F)

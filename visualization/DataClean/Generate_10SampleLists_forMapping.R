library(data.table)
library(magrittr)

dt.meta <- readRDS('./RDS/metadata.rds')
Lst.filenames <- dt.meta$Colnames.miRNA

i = 3


for(i in 1:10){
  
  start = (i-1)*32+1
  end = i*32
  
  Lst.filenames.i <- Lst.filenames[start:end]
  ID.SampleList <- sprintf('SampleList_%d.txt', i)
  dir.SampleList <- paste0('./SampleLists/', ID.SampleList)
  write.table(Lst.filenames.i, dir.SampleList, sep = '\t', quote = F, row.names = F, col.names = F)
  
  print(i)
}

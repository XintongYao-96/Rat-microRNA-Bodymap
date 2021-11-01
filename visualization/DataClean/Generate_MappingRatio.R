rm(list = ls())

library(data.table)
library(magrittr)
dt.exprSet <- readRDS('./RDS/exprMat/exprSet_rawcount_1033miR320Sample.rds')

# Genome stats ------------------------------------------------------------

dir.readcount <- './ReadStats_Log/'
folders.readcount <- dir(dir.readcount)

x <- folders.readcount[1]
dt.ReadStats <- lapply(folders.readcount, function(x){
  
  dir.currFile <- sprintf('%s%s', dir.readcount, x,x)
  dt.in <- fread(dir.currFile, fill = T)
  dt.in$TotalSeq <- lapply(strsplit(as.character(dt.in$V1), '\\:'), '[', 2)
  
  Total.Genome <- dt.in[1,]$TotalSeq
  Genome.mapped <- dt.in[1,2]
  miRNA.mapped <- dt.in[2,2]
  Sample.ID = strsplit(x, split = "\\.")[[1]][1]
  
  dt.stats = data.table(Sample.ID = Sample.ID,
                        Genome.Mapped = Genome.mapped,
                        miRNA.Mapped = miRNA.mapped)

  return(dt.stats)
})%>%rbindlist()

setnames(dt.ReadStats, c('Sample.ID', 'Genome.Mapped', 'miRNA.Mapped'))



# miRNA stats -------------------------------------------------------------

IDs.miRNA <- unique(dt.exprSet$miRNA.ID)
IDs.novel.miRNA <- grep('chr', IDs.miRNA, value = T)

# novel miRNA Reads
dt.exprSet.novel <- dt.exprSet[miRNA.ID%in%IDs.novel.miRNA]
dt.N.miRNA.novel <- dt.exprSet.novel[, .(Novel.miRNA.Mapped = round(sum(Count))), by = .(Sample.ID)]

# dt.exprSet.known <- dt.exprSet[miRNA.ID%in%IDs.known.miRNA]
# dt.N.miRNA.known <- dt.exprSet.known[, .(Known.miRNA.Mapped = round(sum(Count))), by = .(Sample.ID)]

# dt.miRNA.ReadStats <- merge(dt.N.miRNA.known, dt.N.miRNA.novel)

dt.merge.ReadStats <- merge(dt.ReadStats, dt.N.miRNA.novel, by = 'Sample.ID')
dt.merge.ReadStats$Known.miRNA.Mapped <- dt.merge.ReadStats$miRNA.Mapped-(dt.merge.ReadStats$Novel.miRNA.Mapped)

saveRDS(dt.merge.ReadStats, './RDS/ReadStats.rds')

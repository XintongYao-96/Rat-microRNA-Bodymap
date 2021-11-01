rm(list = ls())
library(magrittr)
library(data.table)

dt.meta <- fread('./tables/Table1_metadata_forEdit.csv')
dt.ReadStats <- readRDS('./RDS/ReadStats.rds')
dt.N.miRNA <- readRDS('./RDS/metadata_with_miRNA_Number.rds')

# dt.ReadStats.forMerge <- dt.ReadStats[, .(Sample.ID, Genome.mapped, Genome.unmapped, Total.known.miRNA, Total.novel.miRNA, miRNA.unmapped)]
dt.N.miRNA.forMerge <- dt.N.miRNA[, .(ID.sample, Total.miRNA.Number, Known.miRNA.Number, Novel.miRNA.Number)]

dt.meta.merge <- merge(dt.meta, dt.ReadStats, by.x = 'Colnames.miRNA', by.y = 'Sample.ID')
dt.meta.merge$Genome.Reads.Proportion <- dt.meta.merge$Genome.Mapped/(dt.meta.merge$`Total sequence`)
dt.meta.merge$miRNA.Reads.Proportion <- (dt.meta.merge$Known.miRNA.Mapped+dt.meta.merge$Novel.miRNA.Mapped)/dt.meta.merge$`Total sequence`

dt.meta.merged <- merge(dt.meta.merge, dt.N.miRNA.forMerge, by.x = 'Sample_ID', by.y = 'ID.sample')

write.csv(dt.meta.merged, './tables/Table1_metadata_Edited.csv', row.names = F)


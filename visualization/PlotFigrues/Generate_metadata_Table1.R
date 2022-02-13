rm(list = ls())
library(magrittr)
library(data.table)

dt.meta <- fread('./tables/Table1_metadata_forEdit.csv')
dt.readstats <- readRDS('./RDS_202201/ReadStats.rds')

dt.readstat.dcast <- dcast(dt.readstats, Sample.ID~Type, value.var = 'Count')


dt.meta.merge <- merge(dt.meta, dt.readstat.dcast, by.x = 'Colnames.miRNA', by.y = 'Sample.ID')
# write.csv(dt.meta.merge, './tables/Table1_Metadata.csv', row.names = F)
dt.sra <- fread('./tables/SraRunTable (4).txt')
dt.GSM <- fread('./tables/GSM_sampleID_link.csv')


dt.SRA.GSM <- data.table(GEO_Accession=dt.sra$`GEO_Accession (exp)`,
                         SRA_ID=dt.sra$Run)
dt.SRA.GSM.SampleID <- merge(dt.SRA.GSM, dt.GSM, by.x  = 'GEO_Accession', by.y ='GSM' )

dt.meta.forPublish <- merge(dt.meta.merge,dt.SRA.GSM.SampleID, by = 'Sample_ID' )

write.csv(dt.meta.forPublish, './tables/Table1_metadata_forPublish_220207.csv')



# stats for the outlier brain sample --------------------------------------

dt.meta.merge[Colnames.miRNA=='lane5_Brn_F_100_4_1']


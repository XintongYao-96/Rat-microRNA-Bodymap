source('./scripts/library.R')
source('./scripts/parameters.R')

dt.DEGresult <- readRDS('./RDS/DEG_results/DEG_OrganPerAge_combine.rds')

# find organ enriched miRNAs
cutOff=4
dt.miRNAs.organEnriched.FC4 <- dt.DEGresult[adj.P.Val<0.05][logFC>(log2(cutOff))][,.N,by=.(Organ.B,ID.gene)][N==40]
setnames(dt.miRNAs.organEnriched.FC4,'Organ.B','Organ')

#dt.miRNAs.organEnriched[, .(n= .N), by= 'Organ']

cutOff=3
dt.miRNAs.organEnriched.FC3 <- dt.DEGresult[adj.P.Val<0.05][logFC>(log2(cutOff))][,.N,by=.(Organ.B,ID.gene)][N==40]
setnames(dt.miRNAs.organEnriched.FC3,'Organ.B','Organ')

cutOff=2
dt.miRNAs.organEnriched.FC2 <- dt.DEGresult[adj.P.Val<0.05][logFC>(log2(cutOff))][,.N,by=.(Organ.B,ID.gene)][N==40]
setnames(dt.miRNAs.organEnriched.FC2,'Organ.B','Organ')






# find organ lacked miRNAs
cutOff=2
dt.miRNAs.organLacked <- dt.DEGresult[adj.P.Val<0.05][logFC< -(log2(cutOff))][,.N,by=.(Organ.B,ID.gene)][N==40]
setnames(dt.miRNAs.organLacked,'Organ.B','Organ')


# write tables

write.table(dt.miRNAs.organEnriched, './tables/list_organEnriched_miRNAs.txt', row.names = F, quote = F)
write.table(dt.miRNAs.organLacked, './tables/list_organLacked_miRNAs.txt', row.names = F, quote = F)
saveRDS(dt.miRNAs.organLacked, './RDS/dt.miRNAs.organLacked.rds')
saveRDS(dt.miRNAs.organEnriched, './RDS/dt.miRNAs.organEnriched')



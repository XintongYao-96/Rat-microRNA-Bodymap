dt.DEGresult <- readRDS('./RDS/DEG_results/DEG_OrganPerAge_combine.rds')


lst.cutoff <- c(1,2,4,6,8,10)
dt.miRNA.OrganEnriched.multiCutoff <- sapply(lst.cutoff, function(x){
  
  cutoff <- x
  dt.miRNAs.organEnriched <- dt.DEGresult[adj.P.Val<0.05][logFC>(log2(cutoff))][,.N,by=.(Organ.B,ID.gene)][N==40]
  dt.miRNAs.organEnriched$cutoff <- cutoff
  
  return(list(dt.miRNAs.organEnriched))
  print((log2(cutoff)))
})%>%rbindlist()


dt.miRNA.OrganEnriched.result <- data.table(Organ = dt.miRNA.OrganEnriched.multiCutoff$Organ.B,
                                            miRNA = dt.miRNA.OrganEnriched.multiCutoff$ID.gene, 
                                            cutoff = dt.miRNA.OrganEnriched.multiCutoff$cutoff)

write.csv(dt.miRNA.OrganEnriched.result, './charts/SupplementaryTables/OrganEnriched_miRNAs_multicutoff.csv', row.names = F)

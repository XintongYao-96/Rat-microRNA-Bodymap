source('./scripts/library.R')
source('./scripts/Functions/Function.Call.DEGs.cpm.R')
source('./scripts/parameters.R')
exprMat.rawCount <- readRDS('./RDS/exprMat/exprMat_miRNA_rawFilt.rds')

organs.sexUnrelated <- setdiff(organs,c('Tst','Utr'))


dt.DEGanaly <- rbindlist(lapply(organs.sexUnrelated, function(currOrgan){
    
    print(currOrgan)
    
    # Call DEGs between genders of samples from current organ in each development stage
    dt.DEGs.currOrgan <- rbindlist(lapply(ages,function(currAge){
      
      IDs.sample.currGroup <- c(dt.meta[Organ==currOrgan][Age==currAge][Sex=='F']$ID.sample,
                                dt.meta[Organ==currOrgan][Age==currAge][Sex=='M']$ID.sample)
      IDs.sample.subset <- intersect(colnames(exprMat.rawCount),IDs.sample.currGroup)
      
      exprMat <- exprMat.rawCount[,IDs.sample.subset]
      groupInfo <- data.frame(ID.sample=IDs.sample.subset,
                              Group=as.character(dt.meta[IDs.sample.subset]$Sex),
                              row.names=1)
      
      dt.DEGs.currOrgan.currAge <- Func.CallDEGs(exprMat,groupInfo)
      
      dt.DEGs.currOrgan.currAge$ID.gene <- rownames(dt.DEGs.currOrgan.currAge)
      dt.DEGs.currOrgan.currAge$Age <- currAge
      dt.DEGs.currOrgan.currAge$Organ <- currOrgan
      return(dt.DEGs.currOrgan.currAge)
    }))
    return(dt.DEGs.currOrgan)
  }))

dt.DEGanaly.lowFilt <- dt.DEGanaly[AveExpr>=1]

saveRDS(dt.DEGanaly.lowFilt,'./RDS/DEG_results/dt.DEG.SexPerOrganAge.cpm.rds')

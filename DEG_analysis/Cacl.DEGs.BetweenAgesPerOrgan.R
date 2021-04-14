source('./scripts/library.R')
source('./scripts/parameters.R')
source('./scripts/Functions/Function.Call.DEGs.cpm.R')

dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.rawCount <- readRDS('./RDS/exprMat/exprMat_miRNA_rawFilt.rds')

pairs.ages <- t(combn(as.character(ages),m=2))

dt.DEGanaly <- rbindlist(
  x<-apply(pairs.ages,1,function(pair.currAge){
    age.A <- pair.currAge[1]
    age.B <- pair.currAge[2]
    
    print(paste(age.A,age.B,sep=' vs '))
    
    # Call DEGs between the pair per age
    dt.DEGs.currAgePair <- rbindlist(lapply(organs,function(currOrgan){
      
      IDs.sample.currGroup <- c(dt.meta[Organ==currOrgan][(Age==age.A)]$ID.sample,
                                dt.meta[Organ==currOrgan][(Age==age.B)]$ID.sample)
      IDs.sample.subset <- intersect(colnames(exprMat.rawCount),IDs.sample.currGroup)
      
      exprMat <- exprMat.rawCount[,IDs.sample.subset]
      groupInfo <- data.frame(ID.sample=IDs.sample.subset,
                              Group=as.character(dt.meta[IDs.sample.subset]$Age),
                              row.names=1)
      
      dt.DEGs.currAgePair.currOrgan <- Func.CallDEGs(exprMat,groupInfo)
      
      dt.DEGs.currAgePair.currOrgan$ID.gene <- rownames(dt.DEGs.currAgePair.currOrgan)
      dt.DEGs.currAgePair.currOrgan$Age.A <- min(age.A,age.B)
      dt.DEGs.currAgePair.currOrgan$Age.B <- max(age.A,age.B)
      dt.DEGs.currAgePair.currOrgan$Organ <- currOrgan
      return(dt.DEGs.currAgePair.currOrgan)
    }))
    return(dt.DEGs.currAgePair)
  }))

# Filter DEGs from two low-exprssed groups

dt.DEGanaly.lowFilt <- dt.DEGanaly[AveExpr>=1]

saveRDS(dt.DEGanaly.lowFilt,'./RDS/DEG_results/DEG_AgeperOrgan.rds')



rm(list = ls())
source('./scripts/parameters.R')
source('./scripts/Functions/Function.Call.DEGs.cpm.R')

dt.meta <- readRDS('./RDS/metadata.rds')
exprMat.rawCount <- readRDS('./RDS/exprMat_miRNA_rawFilt.rds')

pairs.organs <- t(combn(as.character(organs),m=2))

dt.DEGanaly <- rbindlist(
  x<-apply(pairs.organs,1,function(pair.currOrgan){
    organ.A <- pair.currOrgan[1]
    organ.B <- pair.currOrgan[2]

    print(paste(organ.A,organ.B,sep=' vs '))
 
    # Call DEGs between the pair per age
    dt.DEGs.currOrganPair <- rbindlist(lapply(ages,function(currAge){
  
        IDs.sample.currGroup <- c(dt.meta[Age==currAge][(Organ==organ.A)]$ID.sample,dt.meta[Age==currAge][(Organ==organ.B)]$ID.sample)
        IDs.sample.subset <- intersect(colnames(exprMat.rawCount),IDs.sample.currGroup)
  
        exprMat <- exprMat.rawCount[,IDs.sample.subset]
        groupInfo <- data.frame(ID.sample=IDs.sample.subset,
                                Group=as.character(dt.meta[IDs.sample.subset]$Organ),
                                row.names=1)

        dt.DEGanaly.currOrganpair.currAge <- Func.CallDEGs(exprMat,groupInfo)
        
        dt.DEGanaly.currOrganpair.currAge$ID.gene <- rownames(dt.DEGanaly.currOrganpair.currAge)
        dt.DEGanaly.currOrganpair.currAge$Age <- currAge
        dt.DEGanaly.currOrganpair.currAge$Organ.A <- organ.A
        dt.DEGanaly.currOrganpair.currAge$Organ.B <- organ.B
        return(dt.DEGanaly.currOrganpair.currAge)
      }))
    return(dt.DEGs.currOrganPair)
}))


# Filter DEGs from two low-exprssed groups

dt.DEGanaly.lowFilt <- dt.DEGanaly[AveExpr>=1]
dt.DEGs <- dt.DEGanaly.lowFilt

# combine th OrganB v.s. A result with OrganB v.s. A result 
dt.result.rev <- data.table(logFC=-dt.DEGs$logFC,
                            AveExpr=dt.DEGs$AveExpr,
                            t=dt.DEGs$t,
                            P.Value=dt.DEGs$P.Value,
                            adj.P.Val=dt.DEGs$adj.P.Val,
                            B=dt.DEGs$B,
                            ID.gene=dt.DEGs$ID.gene,
                            Age=dt.DEGs$Age,
                            Organ.A=dt.DEGs$Organ.B,
                            Organ.B= dt.DEGs$Organ.A)

dt.results <- list(dt.DEGs, dt.result.rev)%>%rbindlist()

saveRDS(dt.DEGanaly.lowFilt,'./RDS/DEG_OrganPerAge_org.rds')
saveRDS(dt.results,'./RDS/DEG_OrganPerAge_combine.rds')



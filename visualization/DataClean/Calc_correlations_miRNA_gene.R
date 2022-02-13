library(data.table)
library(magrittr)

dt.cor.result <- readRDS('./RDS/Corr_allType.rds')
exprMat.miRNA <- readRDS('./RDS/exprMat/exprMat_logCPM_r604c318.rds')
exprMat.mRNA <- readRDS('./RDS/exprMat/exprMat_mRNA_log2FPKM.rds')


IDs.samples <- intersect(colnames(exprMat.miRNA), colnames(exprMat.mRNA))

exprMat.bind <- rbind(exprMat.miRNA[, IDs.samples], exprMat.mRNA[, IDs.samples])

corr.mat <- cor(t(exprMat.bind))
corr.miRNA.mRNA <- corr.mat[rownames(exprMat.miRNA), rownames(exprMat.mRNA)]

dt.melt.corr <- melt(corr.miRNA.mRNA)
dt.melt.corr$pair <- paste0(dt.melt.corr$ID.Gene.A, dt.melt.corr$ID.Gene.B)
# setnames(dt.melt.corr, c('ID.Gene.A', 'ID.Gene.B', 'corr'))
dt.melt.corr <- data.table(dt.melt.corr)


# miRNA- host gene (same)

dt.host.same.pairs <- dt.cor.result[Type == 'miRNA-host gene (same)']
dt.host.same.pairs$pair <- paste0(dt.host.same.pairs$ID.Gene.A,dt.host.same.pairs$ID.Gene.B )
dt.corr.miRNA.hostGene.same <- dt.melt.corr[pair%in%dt.host.same.pairs$pair]

# miRNA-host gene (oppo)
dt.host.oppo.pairs <- dt.cor.result[Type == 'miRNA-host gene (opp)']
dt.host.oppo.pairs$pair <- paste0(dt.host.oppo.pairs$ID.Gene.A,dt.host.oppo.pairs$ID.Gene.B )
dt.corr.miRNA.hostGene.oppo <- dt.melt.corr[pair%in%dt.host.oppo.pairs$pair]

setdiff(dt.host.oppo.pairs$pair , dt.corr.miRNA.hostGene.oppo$pair)
# combine

dt.same <- data.table(ID.Gene.A = dt.corr.miRNA.hostGene.same$ID.Gene.A,
                      ID.Gene.B = dt.corr.miRNA.hostGene.same$ID.Gene.B,
                      corr = dt.corr.miRNA.hostGene.same$corr,
                      Type = 'miRNA-host gene (same)')
dt.oppo <- data.table(ID.Gene.A = dt.corr.miRNA.hostGene.oppo$ID.Gene.A,
                      ID.Gene.B = dt.corr.miRNA.hostGene.oppo$ID.Gene.B,
                      corr = dt.corr.miRNA.hostGene.oppo$corr,
                      Type = 'miRNA-host gene (oppo)')
dt.all <- dt.melt.corr[, .(ID.Gene.A, ID.Gene.B, corr)]
dt.all$Type <- 'miRNA-all gene'

dt.corr.combined <- rbind(dt.same, dt.oppo, dt.all)

saveRDS(dt.corr.combined, './RDS/Corr_miRNA_Gene_combined.rds')

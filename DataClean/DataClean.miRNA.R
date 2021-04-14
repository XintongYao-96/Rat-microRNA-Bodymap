rm(list = ls())
library(data.table)
exprMat.raw <- as.matrix(read.table('./data/RawData/bodymap_miRNA_reads_count.txt'))

# ID transfer -------------------------------------------------------------
dt.meta <- readRDS('RDS/metadata.rds')
setkey(dt.meta,Colnames.miRNA)
colnames(exprMat.raw) <- dt.meta[colnames(exprMat.raw)]$ID.sample


# miRNA filtration --------------------------------------------------------

# Criteria：
# 1）Drop non-chormosomal novel miRNAs
# 2) Drop consistently low-expressed miRNAs

IDs.miRNA.all <- rownames(exprMat.raw)
IDs.miRNA.nonChorm <- IDs.miRNA.all[grep('chrUn|chrM',IDs.miRNA.all)]
IDs.miRNA.lowExpr <- rownames(exprMat.raw[apply(exprMat.raw,1,sum)==0,])
IDs.miRNA.passFilt <- setdiff(IDs.miRNA.all,c(IDs.miRNA.nonChorm,IDs.miRNA.lowExpr))

exprMat.mirFilt <- exprMat.raw[IDs.miRNA.passFilt,]

# Write fitler logs
dt.miRNA.stage <- data.table(
  ID.miRNA = c(IDs.miRNA.nonChorm,IDs.miRNA.lowExpr,IDs.miRNA.passFilt),
  Stage = c(rep('FILTERED',length(IDs.miRNA.nonChorm)+length(IDs.miRNA.lowExpr)),
             rep('PASS',length(IDs.miRNA.passFilt))
             ),
  Comment = c(rep('non-choromsomal novel miRNAs',length(IDs.miRNA.nonChorm)),
               rep('no expression in all samples',length(IDs.miRNA.lowExpr)),
               rep('',length(IDs.miRNA.passFilt))
               )
)

write.csv(dt.miRNA.stage,'./tables/miRNA filtration.csv',row.names = FALSE, quote=FALSE)


# Drop QC-failed samples --------------------------------------------------

# Criteria:
# 1) drop samples failed to cluster with others in group

IDs.sample.all <- colnames(exprMat.raw)
IDs.sample.failedCluster <- c('Spl_F_104_4','Brn_M_006_3')
IDs.sample.passFilt <- setdiff(IDs.sample.all,IDs.sample.failedCluster)

exprMat.lowFilt.sampleFilt <- exprMat.mirFilt[,IDs.sample.passFilt]

write.table(exprMat.lowFilt.sampleFilt,
            sprintf('./tables/ExprMat_miRNA_rawFilt_%dr%dc.txt',nrow(exprMat.lowFilt.sampleFilt),ncol(exprMat.lowFilt.sampleFilt)),
            quote=FALSE)

saveRDS(exprMat.lowFilt.sampleFilt,'./RDS/exprMat_miRNA_rawFilt.rds')



# Normalization -----------------------------------------------------------

library(edgeR)

setkey(dt.meta,ID.sample)
dge <- DGEList(exprMat.lowFilt.sampleFilt,group = as.numeric(dt.meta[colnames(exprMat.lowFilt.sampleFilt)]$Group))
dge <- calcNormFactors(dge, method='TMM')
exprMat.cpm <- cpm(dge, normalized.lib.sizes = FALSE, log=FALSE, prior.count = 1)
exprMat.logCPM <- cpm(dge, normalized.lib.sizes = TRUE, log=TRUE, prior.count = 1)


# Save DGE
saveRDS(dge,'./RDS/dge.rds')


# Save normalized matrix
write.table(exprMat.logCPM,
            sprintf('./tables/ExprMat_miRNA_TMMlog2CPM_%dr%dc.txt',nrow(exprMat.lowFilt.sampleFilt),ncol(exprMat.lowFilt.sampleFilt)),
            quote=FALSE)

saveRDS(exprMat.logCPM,'./RDS/exprMat_miRNA_log2CPM.rds')
saveRDS(exprMat.cpm, './RDS/exprMat/exprMat_miRNA_CPM.rds')

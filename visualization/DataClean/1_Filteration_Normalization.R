rm(list = ls())
library(data.table)
library(magrittr)
library(edgeR)

exprMat.raw <- readRDS('./RDS/exprMat/exprMat_rawcount_r745c320.rds')


# miRNA filteration -------------------------------------------------------
# drop constently low-expressed miRNAs (low-expressed: count<20 in 320 samples)

IDs.miRNA.all <- rownames(exprMat.raw)

IDs.miRNA.lowExpr <- rownames(exprMat.raw[apply(exprMat.raw,1,sum)<=320,])
IDs.miRNA.passFilt <- setdiff(IDs.miRNA.all,IDs.miRNA.lowExpr)
length(IDs.miRNA.passFilt)
exprMat.miRFlit <- exprMat.raw[IDs.miRNA.passFilt,]

# count.miRNA.merge <- rowSums(exprMat.miRFlit)

# grep('chr', IDs.miRNA.passFilt)%>%length()


# ID transfer -------------------------------------------------------------

dt.meta <- readRDS('RDS/metadata.rds')


# Sample filteration ------------------------------------------------------
# Drop 2 samples fail to cluster with the samples from given organs
IDs.sample.all <- colnames(exprMat.miRFlit)
IDs.sample.failedCluster <- c("Brn_M_006_3", "Spl_F_104_4")
IDs.sample.passFilt <- setdiff(IDs.sample.all,IDs.sample.failedCluster)

exprMat.miRFilt.sampleFilt <- exprMat.miRFlit[,IDs.sample.passFilt]


# Normalization -----------------------------------------------------------

setkey(dt.meta,ID.sample)
dge <- DGEList(exprMat.miRFilt.sampleFilt,group = as.numeric(dt.meta[colnames(exprMat.miRFilt.sampleFilt)]$Group))
dge <- calcNormFactors(dge, method='TMM')
exprMat.cpm <- cpm(dge, normalized.lib.sizes = FALSE, log=FALSE, prior.count = 1)
exprMat.logCPM <- cpm(dge, normalized.lib.sizes = TRUE, log=TRUE, prior.count = 1)



# save --------------------------------------------------------------------
write.table(exprMat.raw, './tables/GSE172269_RatBodyMap_miRNA_RawCount.txt', quote = F, sep = '\t')
saveRDS(exprMat.miRFilt.sampleFilt, './RDS/exprMat/exprMat_rawcount_r604c318.rds')
saveRDS(exprMat.cpm, './RDS/exprMat/exprMat_CPM_r604c318.rds')
saveRDS(exprMat.logCPM, './RDS/exprMat/exprMat_logCPM_r604c318.rds')

# expmat.geo <- fread('./tables/GSE172269_RatBodyMap_miRNA_RawCount.txt')
# # 
# # Test HCA ---------------------------------------------------------------------
# # Perform 320 samples (no-filter) HCA
# 
source("./scripts/PlotFigures/parameters.R")
library(pheatmap)

exprMat.forHCA <- exprMat.logCPM

annot_col <- data.frame(dt.meta[,.(ID.sample,Sex,Age,Organ)],row.names=1)
annot_col$Sex <- factor(annot_col$Sex, levels = c('M', 'F'))

color.cluster.type <- c('#d9ecf2','#f56a79')
names(color.cluster.type) <- c('none', 'cluster1')

annot_color <- list(Organ=colors.organ,
                    Age=colors.age,
                    type=color.cluster.type,
                    Sex=color.sex )

pheatmap(exprMat.forHCA,
         show_rownames = FALSE, show_colnames = FALSE,
         treeheight_row = 40, treeheight_col = 60,
         clustering_method = "ward.D",scale="row",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[9:11]))(29),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[3:9]))(59),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:3]))(29)),
         breaks = c(seq(-17,-3,0.5),
                    seq(-2.9,2.9,0.1),
                    seq(3,17,0.5)),
         annotation_col = annot_col, annotation_colors = annot_color,
         annotation_names_col=FALSE)


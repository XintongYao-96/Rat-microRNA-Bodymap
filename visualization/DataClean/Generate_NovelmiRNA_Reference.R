#######################################
#       Novel miRNA Exploration       #
#                                     #
#            Xintong Yao              #
#              2022-01                #
#                                     #
# Get fasta reference                 #
# Combine novel exprmat               #
# Report novel miRNA  seq (Table 2）  #
#                                     #
#######################################

rm(list = ls())
library(data.table)
library(magrittr)
library(RColorBrewer)
library(pheatmap)


dt.miRDeep.result <- fread('./tables/miRDeepTable_NovelMiRNA_merge320.csv')
dt.meta <- readRDS('./RDS/metadata.rds')

# Primary Filtering Criteria: score >5 and p-value<0.05 
dt.pass.novel.miRNA <- dt.miRDeep.result[`miRDeep2 score`>=4][`significant randfold p-value`=='yes'][`rfam alert`!= 'rRNA']
IDs.novel.miRNA <- paste0('rno_miR_', dt.pass.novel.miRNA$`provisional id`)
IDs.novel.premiRNA <- paste0('rno_mir_', dt.pass.novel.miRNA$`provisional id`)

dt.novel.info <- data.table(miRNA_ID = IDs.novel.miRNA,
                            miRDeep2_score = dt.pass.novel.miRNA$`miRDeep2 score`,
                            mature_seq = toupper(dt.pass.novel.miRNA$`consensus mature sequence`),
                            star_seq = toupper(dt.pass.novel.miRNA$`consensus star sequence`),
                            pre_seq = toupper(dt.pass.novel.miRNA$`consensus precursor sequence`),
                            position = toupper(dt.pass.novel.miRNA$`precursor coordinate`)
)


# Get fasta reference -----------------------------------------------------
# write fasta
SEQs.novel.miRNA <- gsub('U', 'T', toupper(dt.pass.novel.miRNA$`consensus mature sequence`))
SEQs.novel.pre.miRNA <- gsub('U', 'T', toupper(dt.pass.novel.miRNA$`consensus precursor sequence`))


N.rows = 2*nrow(dt.novel.info)
Index.odd <- seq(1, N.rows-1, 2)
Index.even <- seq(2, N.rows, 2)

# mature miRNA
Lst.reference <- rep(0, N.rows)
Lst.reference[Index.odd]<- paste0('>', IDs.novel.miRNA)
Lst.reference[Index.even]<- SEQs.novel.miRNA
head(Lst.reference)
write.table(Lst.reference, './tables/rno_mature_novel.fa',sep = '\t', quote = F, row.names = F, col.names = F)

# precursor miRNA
Lst.premiRNA.ref <- rep(0, N.rows)
Lst.premiRNA.ref[Index.odd]<- paste0('>', IDs.novel.premiRNA)
Lst.premiRNA.ref[Index.even]<- SEQs.novel.pre.miRNA
write.table(Lst.premiRNA.ref, './tables/rno_precursor_novel.fa',sep = '\t', quote = F, row.names = F, col.names = F)


# Report --------------------------------------------------------------

novel.exprmat <- readRDS('./RDS/exprMat/exprMat_novel_count_r270c320.rds')


# Annotate piRNA alert
dt.piRBase.map <- fread('./tables/novel_piRBaseAligned.txt')
dt.piRBase.info <- data.table(miRNA_ID = dt.piRBase.map$V1,
                              piRNA_alert = dt.piRBase.map$V3)

dt.novel.info.piRmerge <- merge(dt.novel.info, dt.piRBase.info, by = 'miRNA_ID',
                                all.x = T)


# Organ-specific novel miRNA
dt.DEG <- readRDS('./RDS/Novel_DEG_OrganPerAge_combine.rds')

cutOff=1.5
dt.miRNAs.organEnriched <- dt.DEG[adj.P.Val<0.05][logFC>(log2(cutOff))][,.N,by=.(Organ.B,ID.gene)][N==40]

dt.organ.enriched.novel <- data.table(miRNA_ID = dt.miRNAs.organEnriched$ID.gene,
                                      Organ_enriched = dt.miRNAs.organEnriched$Organ)


# Get the novel miRNA metadata table
dt.organ.enriched.piRNA.annot <- merge(dt.novel.info.piRmerge, dt.organ.enriched.novel, by = 'miRNA_ID', all = F)


# miRDeep score filtering: score >=10
# Filter: remove repeat sequence
lst.repeat.seqs <- dt.organ.enriched.piRNA.annot$mature_seq[table(dt.organ.enriched.piRNA.annot$mature_seq)>1]
dt.novel.annot <- dt.organ.enriched.piRNA.annot[!mature_seq%in%lst.repeat.seqs][miRDeep2_score>=10]

# merge exprmat

dt.exprmat <- data.table(miRNA_ID= rownames(novel.exprmat),
                         novel.exprmat[, dt.meta$ID.sample])

dt.novel.miRNA <- merge(dt.novel.annot, dt.exprmat, by = 'miRNA_ID')

write.csv(dt.novel.miRNA, './Table2_Meta_Novel_miRNA_Information_2022.csv')

#check heatmap
source("./scripts/PlotFigures/parameters.R")
dt.meta <- readRDS('./RDS/metadata.rds')


annot_col <- data.frame(dt.meta[,.(ID.sample,Sex,Age,Organ)],row.names=1)
annot_color <- list(Sex = color.sex,
                    Age=colors.age,
                    Organ=colors.organ)

pheatmap(novel.exprmat[dt.novel.miRNA$miRNA_ID,dt.meta$ID.sample],
         show_rownames = FALSE, show_colnames = FALSE,
         treeheight_col = 60,
         cluster_cols = F,
         clustering_method = "ward.D",
         scale="row",
         color = c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[9:11]))(21),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[3:9]))(39),
                   colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[1:3]))(41)),
         breaks = c(seq(0,2,0.1),
                    seq(2.1,5.9,0.1),
                    seq(6,10,0.1)),
         annotation_col = annot_col, annotation_colors = annot_color,
         annotation_names_col=FALSE,
         filename = './charts/Response_Heatmap_novelmiRNA.pdf', width = 12, height = 8)



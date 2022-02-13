# Expand miRNA Reference 
# 5 nt Upstream and downstream 
# 2022-01-13


# Data clean --------------------------------------------------------------

library(data.table)
library(magrittr)

dt.pos.miRNA <- fread('./reference/miRNA_expand/mature_align2_hairpin.txt')
dt.seq.hairpin <- read.table('./reference/miRNA_expand/rno_hairpin_ref.txt')


# Generate hairpin ID and sequence table

Lst.seq.hairpin <- dt.seq.hairpin$V1

N.rows = length(Lst.seq.hairpin)
Index.odd <- seq(1, N.rows-1, 2)
Index.even <- seq(2, N.rows, 2)
Lst.ID.mir <- Lst.seq.hairpin[Index.odd]%>%substr(., 2, nchar(.))
Lst.seq.mir <- Lst.seq.hairpin[Index.even]
dt.seq.mir <- data.table(ID.mir = Lst.ID.mir, 
                         Seq.mir = Lst.seq.mir)


# Generate miRNA position info for combine


dt.miRNA.pos.forMerge <- data.table(miRNA.ID = dt.pos.miRNA$V1,
                                    ID.mir = dt.pos.miRNA$V3,
                                    Pos = dt.pos.miRNA$V4,
                                    Length.miR = as.numeric(substr(dt.pos.miRNA$V6, 1,2)),
                                    Seq.miR = dt.pos.miRNA$V10)

dt.miRNA.pos.info <- merge(dt.miRNA.pos.forMerge, dt.seq.mir,
                           by = 'ID.mir')


# Expand ------------------------------------------------------------------

# upstream 5 nt
dt.miRNA.pos.info$Expand.start <- ifelse(dt.miRNA.pos.info$Pos>=6, dt.miRNA.pos.info$Pos-5, 1)
# downstream 5 nt
dt.miRNA.pos.info$Expand.end <- 4+dt.miRNA.pos.info$Pos+dt.miRNA.pos.info$Length.miR
# expand sequence
dt.miRNA.pos.info$Seq.miR.expand <- substr(dt.miRNA.pos.info$Seq.mir, dt.miRNA.pos.info$Expand.start, dt.miRNA.pos.info$Expand.end)


# write table -------------------------------------------------------------

Lst.miR.expand <- seq(1, 2*nrow(dt.miRNA.pos.info))
Index.odd <- seq(1, length(Lst.miR.expand)-1, 2) 
Index.even <- seq(2, length(Lst.miR.expand), 2) 

Lst.miR.expand[Index.odd] <- paste0('>', dt.miRNA.pos.info$miRNA.ID)
Lst.miR.expand[Index.even] <- dt.miRNA.pos.info$Seq.miR.expand

write.table(Lst.miR.expand, './reference/miRNA_expand/rno_mature_expand_ref.txt',
            col.names = F, row.names = F, quote = F)

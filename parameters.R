# Dependent libraries
library(data.table)
library(RColorBrewer)
library(ggplot2)

dt.meta <- readRDS('./RDS/metadata.rds')

# Set colors
organs <- unique(dt.meta$Organ)
colors.organ <- colorRampPalette(brewer.pal(9,'Set1'))(11)
names(colors.organ) <- organs

colors.organ.pastel <- colorRampPalette(brewer.pal(9,'Pastel1'))(11)
names(colors.organ.pastel) <- organs
colors.system <- c(`nervous system` = as.character(colors.organ.pastel['Brn']),
                   `urogenital system` = as.character(colors.organ.pastel['Tst']),
                   `immune system` = as.character(colors.organ.pastel['Thm']),
                   `respiratory system` = as.character(colors.organ.pastel['Lng']),
                   `digestive system` = as.character(colors.organ.pastel['Lvr']),
                   `blood circulation system` = as.character(colors.organ.pastel['Hrt']),
                   `endocrine system` = as.character(colors.organ.pastel['Adr']),
                   `motor system` = as.character(colors.organ['Utr'])
                   )


ages <- unique(dt.meta$Age)
colors.age <-c('#bbe1fa','#3282b8','#0f4c75','#1b262c')
names(colors.age) <- ages

sexs <- unique(dt.meta$Sex)
color.sex <- c('#e36387', '#a6dcef')
names(color.sex) <- c('F', 'M')
colors.regulateType <- brewer.pal(3,'Set2')[c(2,1)]
names(colors.regulateType) <- c('Up','Down')

stage.PVCA <- c('Age', 'Sex', 'Age:Sex', 'Residual')%>%rev()
colors.PVCA.stage <- rev(c('#07689f','#a2d5f2','#fafafa','#ff7e67'))
names(colors.PVCA.stage) <- stage.PVCA

systems <- c('Nervous System','Reproductive System','Digestive System',      
             'Immune System','Urinary System','Respiratory System',    
             'Musculoskeletal System', 'Endocrine System','Cardiovascular System')
colors.system <- c('#596A98', "#E485B7", "#FF7F00", "#C9992C" , "#6B886D","#AC5782",
                   "#FFE528", "#E41A1C","#449B75"  )
names(colors.system) <- systems

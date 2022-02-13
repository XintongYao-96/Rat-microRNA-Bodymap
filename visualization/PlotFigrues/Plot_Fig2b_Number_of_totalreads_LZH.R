setwd("~/Downloads/Rat Bodymap/")
library(data.table)
library(magrittr)
library(ggplot2)
# Dependent libraries
library(data.table)
library(RColorBrewer)
library(ggplot2)

# Set colors
dt.meta <- readRDS('./data/Rat_SciData/data/metadata.rds')
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

colors.BioRep <- c('#07689f','#a2d5f2','#fafafa','#ff7e67')
names(colors.BioRep) <- c('1','2','3','4')


stage.PVCA <- c('Age', 'Sex', 'Age:Sex', 'Residual')%>%rev()
colors.PVCA.stage <- rev(c('#07689f','#a2d5f2','#fafafa','#ff7e67'))
names(colors.PVCA.stage) <- stage.PVCA

systems <- c('Nervous System','Reproductive System','Digestive System',      
             'Immune System','Urinary System','Respiratory System',    
             'Musculoskeletal System', 'Endocrine System','Cardiovascular System')
colors.system <- c('#596A98', "#E485B7", "#FF7F00", "#C9992C" , "#6B886D","#AC5782",
                   "#FFE528", "#E41A1C","#449B75"  )
names(colors.system) <- systems


dt.meta <- readRDS('./data/Table_total_read.rds')
colnames(dt.meta)
dt.meta <- as.data.table(dt.meta)
colnames(dt.meta)[20] <- "reads_count"
dt.meta.total_sequences <- dt.meta[, .(Sample_ID, miRNA_count,reads_count)]

dtforPlot <- melt(dt.meta.total_sequences, id.vars = 'Sample_ID',
                  variable.name = 'type_count',
                  value.name = 'Number')
dtforPlot$type_count <- factor(dtforPlot$type_count, levels = c( 'reads_count','miRNA_count'))

color.type_count <- c('#AED3EF', '#2D679B')
dt.meta <- readRDS('./data/Rat_SciData/data/metadata.rds')

dtforPlot <- merge(dt.meta,dtforPlot,by.x="ID.sample",by.y = "Sample_ID")
dtforPlot$order <- 0
list <- unique(dtforPlot$Organ)
for (i in 1:length(list)) {dtforPlot$order[which(dtforPlot$Organ==list[i])] <- LETTERS[i]}

my_facet <- facet_grid(~ order, scales = 'free_x', space = 'free_x')
# levels(dt.forPlot_1$batch) <- dt.forPlot_sum$batch
# levels(dt.forPlot_1$libraryPrep.x) 
# levels(dt.forPlot_1$sum) <- unique(dt.forPlot_1$sum)[order(unique(dt.forPlot_1$sum),decreasing = T)]
dtforPlot$BioRep <- as.factor(dtforPlot$BioRep)

p.main <- ggplot(dtforPlot, aes(x=ID.sample, y=Number,fill=type_count))+ 
  geom_bar(stat="identity",width=0.5,position='dodge')+
  facet_grid(~ order, scales = 'free_x', space = 'free_x')+
  # my_facet+
  theme_few()+
  theme(panel.border = element_blank(),
        strip.text.x = element_text(size=0),
        panel.spacing.x = unit(0.2, "mm")
        # strip.background = element_rect(color='black')
  )+
  labs(x='')+
  theme(axis.text.x = element_blank(),
        #axis.text.x =  element_text(angle = 90),
        axis.ticks = element_blank())+
  scale_fill_manual(values=color.type_count);p.main

p.bar.organ <- ggplot(dtforPlot, aes(x=ID.sample,y=1))+
  geom_bar(stat='identity', aes(fill=Organ),width=1)+
  theme_void()+
  my_facet+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=0)
  )+
  scale_fill_manual(values = colors.organ); p.bar.organ

p.bar.sex <- ggplot(dtforPlot, aes(x=ID.sample,y=1))+
  geom_bar(stat='identity', aes(fill=Sex),width=1)+
  theme_void()+
  my_facet+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=0)
  )+
  scale_fill_manual(values = color.sex); p.bar.sex

p.bar.age <- ggplot(dtforPlot, aes(x=ID.sample,y=1))+
  geom_bar(stat='identity', aes(fill=Age),width=1)+
  my_facet+
  theme_void()+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=0)
  )+
  scale_fill_manual(values = colors.age); p.bar.age 


p.bar.BioRep <- ggplot(dtforPlot, aes(x=ID.sample,y=1))+
  geom_bar(stat='identity', aes(fill=BioRep),width=1)+
  my_facet+
  theme_void()+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=0)
  )+
  scale_fill_manual(values = colors.BioRep); p.bar.BioRep


# plog main figures
plot <- plot_grid(p.main+theme(legend.position = 'none'),
                  p.bar.organ+theme(legend.position = 'none'),
                  p.bar.sex+theme(legend.position = 'none'),
                  p.bar.age+theme(legend.position = 'none'),
                  p.bar.BioRep+theme(legend.position = 'none'),
                  align = "v", ncol = 1, 
                  axis = "tb", 
                  rel_heights = c(15,0.8,0.8,0.8,0.8,0.8)); plot
legend <- plot_grid(get_legend(p.main), get_legend(p.bar.organ), get_legend(p.bar.sex), ncol = 3)                    

all_plot <- plot_grid(plot, legend, nrow = 1, rel_widths = c(10,6)); all_plot

ggsave('./plot/Fig2b_N.pdf',all_plot, width = 14,height = 4)

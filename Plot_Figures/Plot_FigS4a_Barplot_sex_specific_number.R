source('./scripts/parameters.R')

dt.DEGs <- readRDS('./RDS/DEG_results/dt.DEG.SexPerOrganAge.cpm.rds')



# Barplot of numbers of DEG in each organ in each devolopment stage -------


dt.DEGs[,Type:=ifelse(adj.P.Val>=0.05,'Not specific',
                      ifelse(logFC>1, 'Male specific',
                             ifelse(logFC<(-1), 'Female specific',
                                    'Not specific')))]


ct.DEGs <- dt.DEGs[,.N,by=.(Organ,Age,Type)]

dt.forPlot <- ct.DEGs[Type!='Not specific']
dt.forPlot[,Value:=ifelse(Type=='Male specific',N,-N)]

dt.forPlot$Age <- factor(dt.forPlot$Age)
dt.forPlot$Type <- factor(dt.forPlot$Type, c('Male specific','Female specific') )


ggplot(dt.forPlot,aes(x=Age,y=Value))+
  geom_bar(stat = 'identity',pos='stack',color='gray80',
           aes(fill=Type))+
  # geom_point()+
  theme_bw()+
  labs(x='Age',y='Number of sex-specific miRNAs',
       fill='')+
  theme(axis.text=element_text(size=10,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.position = c(0.9,0.3),
        legend.text = element_text(size=12,face='bold'),
        strip.background = element_rect(color='black',fill='white'),
        strip.text = element_text(size=14,face='bold')
  )+
  scale_y_continuous(limits = c(-50, 50))+
  scale_fill_manual(values=brewer.pal(3,'Set2')[c(2,1)])+
  facet_wrap(~Organ,nrow = 2)
ggsave('./charts/Fig.5a.Numbers of sex-specific genes in each organ per age.pdf',width=10,height=5) 



# volcano plot ------------------------------------------------------------


colors.DEG <- c("#7DC0A6",'#ED6F51','#CDC9C4')

dt.DEGs$Type <- factor(dt.DEGs$Type)
ggplot(dt.DEGs, aes(x=logFC, y = -log10(P.Value)))+
  geom_point(aes(color = Type, shape = dt.DEGs$Organ))+
  scale_shape_manual(values = seq(1,9))+
  facet_wrap(~Age, nrow = 1)+
  scale_color_manual(values = colors.DEG)+
  theme_bw()+
  scale_x_continuous(limits = c(-5,5))+
  labs(x = 'log2 (fold-change)', y= '-lg (P-value )')+
  theme(axis.text=element_text(size=18,face='bold',color='gray40'),
        axis.title = element_text(size=18,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold'),
        legend.position = 'top',
        strip.background = element_rect(color='black',fill='white'),
        strip.text = element_text(size=24,face='bold')
  )
ggsave('./charts/Fig.S4b.volcano.age.pdf', width = 20, height = 7)
  
  



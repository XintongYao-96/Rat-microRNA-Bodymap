source('./scripts/library.R')
source('./scripts/parameters.R')

dt.DEGresult <- readRDS('./RDS/DEG_results/DEG_OrganPerAge_combine.rds')

dt.DEGs.over <- dt.DEGresult[P.Value<0.05][logFC>1]
dt.DEGs.under <- dt.DEGresult[P.Value<0.05][logFC< -1]

# generate organ and age pairs
dt.organ.pair <- data.table(organ.A = dt.DEGresult$Organ.A,
                            organ.B = dt.DEGresult$Organ.B)%>%unique()

dt.age.pair <- combn(unique(dt.DEGresult$Age),2)%>%t()%>%as.data.table(.)%>%setnames(., c('age.A', 'age.B'))





# calculate Jacard Index --------------------------------------------------

# B v.s. A over expression
dt.Jacard.over <- split(dt.organ.pair, by=c('organ.A','organ.B'))%>%lapply(., function(currOrgan){
  
  #split by each organ pair
  organ.A = currOrgan$organ.A
  organ.B = currOrgan$organ.B
  dt.DEGs.over.currOrgan <- dt.DEGs.over[Organ.A==organ.A][Organ.B==organ.B]
  
  #split each age pair
  dt.Jacard.currOrganpair <- split(dt.age.pair, by=c('age.A', 'age.B'))%>%lapply(., function(currAge){
    
    age.A=currAge$age.A
    age.B=currAge$age.B
    
    lst.ageA.DEGs <- dt.DEGs.over.currOrgan[ Age == age.A]$ID.gene%>%unique()
    lst.ageB.DEGs <- dt.DEGs.over.currOrgan[ Age == age.B]$ID.gene%>%unique()
    length.int <- intersect(lst.ageA.DEGs, lst.ageB.DEGs)%>%length()
    length.uni <- union(lst.ageA.DEGs, lst.ageB.DEGs)%>%length()
    JI <- length.int/length.uni
    
    return(list(
      organ.A=organ.A,
      organ.B=organ.B,
      age.A=age.A,
      age.B=age.B,
      JacardIndex=JI))
  })%>%rbindlist()
  
  return(dt.Jacard.currOrganpair)
  
})%>%rbindlist()



# B v.s. A under expression

dt.Jacard.under <- split(dt.organ.pair, by=c('organ.A','organ.B'))%>%lapply(., function(currOrgan){
  
  #split by each organ pair
  organ.A = currOrgan$organ.A
  organ.B = currOrgan$organ.B
  dt.DEGs.over.currOrgan <- dt.DEGs.under[Organ.A==organ.A][Organ.B==organ.B]
  
  #split each age pair
  dt.Jacard.currOrganpair <- split(dt.age.pair, by=c('age.A', 'age.B'))%>%lapply(., function(currAge){
    
    age.A=currAge$age.A
    age.B=currAge$age.B
    
    lst.ageA.DEGs <- dt.DEGs.over.currOrgan[ Age == age.A]$ID.gene%>%unique()
    lst.ageB.DEGs <- dt.DEGs.over.currOrgan[ Age == age.B]$ID.gene%>%unique()
    length.int <- intersect(lst.ageA.DEGs, lst.ageB.DEGs)%>%length()
    length.uni <- union(lst.ageA.DEGs, lst.ageB.DEGs)%>%length()
    JI <- length.int/length.uni
    
    return(list(
      organ.A=organ.A,
      organ.B=organ.B,
      age.A=age.A,
      age.B=age.B,
      JacardIndex=JI))
  })%>%rbindlist()
  
  return(dt.Jacard.currOrganpair)
  
})%>%rbindlist()


# combine
dt.Jacard.over$DEGtype <- 'over'
dt.Jacard.under$DEGtype <- 'under'
dt.Jacard <- rbind(dt.Jacard.over, dt.Jacard.under)


# dt.Jacard.over.summ <- dt.Jacard.over[, .(n = median(JacardIndex)), by='organ.B']
# dt.Jacard.over.summ <- setorder(dt.Jacard.over.summ, -n)
# level.organ <- dt.Jacard.over.summ$organ.B
level.organ <- c("Brn", "Adr", "Spl", "Lng", "Thm", "Kdn", "Lvr", "Hrt", "Utr",  "Msc", "Tst")



# plot --------------------------------------------------------------------
dt.Jacard$organ.B <- factor(dt.Jacard$organ.B, levels = level.organ)


# total / side by side boxplot 
ggplot(dt.Jacard , aes(x=organ.B, y=JacardIndex,fill=DEGtype))+
  geom_boxplot(aes(fill=organ.B))+
  scale_fill_manual(values = colors.organ)+
  theme_bw()+
  labs(x='Organ',
       y='Jacard Index',
       title = 'Reproducibility of under-expression DEGs')+
  theme(legend.position = 'none',
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5),
        axis.text = element_text(size=12,face='bold'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold'))



# over 
dt.Jacard.over$organ.B <- factor(dt.Jacard.over$organ.B, levels = level.organ)
ggplot(dt.Jacard.over , aes(x=organ.B, y=JacardIndex))+
  geom_boxplot(aes(fill=organ.B))+
  scale_fill_manual(values = colors.organ)+
  theme_bw()+
  labs(x='Organ',
       y='Jaccard Index',
       title = 'Reproducibility of over-expression miRNAs')+
  theme(legend.position = 'none',
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5),
        axis.text = element_text(size=12,face='bold'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold'))
ggsave('./charts/Fig2b.JI.over.pdf', width = 8, height = 3)



# under
ggplot(dt.Jacard.under , aes(x=organ.B, y=JacardIndex))+
  geom_boxplot(aes(fill=organ.B))+
  scale_fill_manual(values = colors.organ)+
  theme_bw()+
  labs(x='Organ',
       y='Jaccard Index',
       title = 'Reproducibility of under-expression miRNAs')+
  theme(legend.position = 'none',
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5),
        axis.text = element_text(size=12,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=12,face='bold'),
        legend.text = element_text(size=12,face='bold',color='gray40'))
ggsave('./charts/Fig2b.JI.under.pdf', width = 6, height = 3)


# numbers for text --------------------------------------------------------

dt.Jacard.over$organ.B <- as.factor(dt.Jacard.over$organ.B)
number.DEG <- dt.Jacard.over[, mean(JacardIndex), by = organ.B]

dt.Jacard.over[organ.B=='Spl'][JacardIndex<0.5]


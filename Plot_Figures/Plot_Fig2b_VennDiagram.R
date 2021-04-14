library(data.table)
library(magrittr)
library(VennDiagram)

dt.DEGresult <- readRDS('./RDS/DEG_results/DEG_OrganPerAge_combine.rds')

dt.Brn.OverExpressed.w2 <- dt.DEGresult[Organ.B=='Brn'][Age==2][P.Value<0.05][logFC>1]$ID.gene%>%unique() 
dt.Brn.OverExpressed.w6 <- dt.DEGresult[Organ.B=='Brn'][Age==6][P.Value<0.05][logFC>1]$ID.gene%>%unique() 
dt.Brn.OverExpressed.w21 <- dt.DEGresult[Organ.B=='Brn'][Age==21][P.Value<0.05][logFC>1]$ID.gene%>%unique() 
dt.Brn.OverExpressed.w104 <- dt.DEGresult[Organ.B=='Brn'][Age==104][P.Value<0.05][logFC>1]$ID.gene%>%unique() 

lst.Brn.overExpressed.total <- c(dt.Brn.OverExpressed.w2, 
                                     dt.Brn.OverExpressed.w6,
                                     dt.Brn.OverExpressed.w21,
                                     dt.Brn.OverExpressed.w104)%>%unique()

venn.diagram( x = list(dt.Brn.OverExpressed.w2,
                       dt.Brn.OverExpressed.w104,
                       dt.Brn.OverExpressed.w6,
                       dt.Brn.OverExpressed.w21
                       ),
              category.names = c("w2" , "w104", 'w6', 'w21'),
              fill = c('cornflowerblue', 'green', 'yellow', 'darkorchid1'),
              col = 'transparent',
             alpha = c(0.5,0.5,0.5,0.5),
             filename = './charts/Fig2b.Venn_brain.png'
             )


dt.Tst.OverExpressed.w2 <- dt.DEGresult[Organ.B=='Tst'][Age==2][P.Value<0.05][logFC>1]$ID.gene%>%unique() 
dt.Tst.OverExpressed.w6 <- dt.DEGresult[Organ.B=='Tst'][Age==6][P.Value<0.05][logFC>1]$ID.gene%>%unique() 
dt.Tst.OverExpressed.w21 <- dt.DEGresult[Organ.B=='Tst'][Age==21][P.Value<0.05][logFC>1]$ID.gene%>%unique() 
dt.Tst.OverExpressed.w104 <- dt.DEGresult[Organ.B=='Tst'][Age==104][P.Value<0.05][logFC>1]$ID.gene%>%unique() 
venn.diagram( x = list(dt.Tst.OverExpressed.w2,
                       dt.Tst.OverExpressed.w104,
                       dt.Tst.OverExpressed.w6,
                       dt.Tst.OverExpressed.w21),
              category.names = c("w2" , "w104", 'w6', 'w21'),
              fill = c('cornflowerblue', 'green', 'yellow', 'darkorchid1'),
              col = 'transparent',
              alpha = c(0.5,0.5,0.5,0.5),
              filename = './charts/venn_testis.png'
)


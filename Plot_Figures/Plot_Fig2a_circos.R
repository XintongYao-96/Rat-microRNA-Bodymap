rm(list = ls())
source('./scripts/library.R')
source('./scripts/parameters.R')
library(circlize)
library(data.table)



# count DEG for each organ pairs and each age -----------------------------

dt.DEGs <- readRDS('./RDS/DEG_results/DEG_OrganPerAge_org.rds')

# calculate number of DEGs

dt.DEGs.calc <- dt.DEGs[P.Value<0.05][logFC>1 | logFC < -1]
length(unique(dt.DEGs.calc$ID.gene))

ct.DEGs.Up.forw <- dt.DEGs[P.Value<0.05][logFC>1][,.(N=-.N),by=.(Organ.A,Organ.B,Age)]
ct.DEGs.Down.forw <- dt.DEGs[P.Value<0.05][logFC<(-1)][,.N,by=.(Organ.A,Organ.B,Age)]
ct.DEGs.forw <- rbindlist(list(ct.DEGs.Up.forw,ct.DEGs.Down.forw))

ct.DEGs.rev <- data.table(
  Organ.A = ct.DEGs.forw$Organ.B,
  Organ.B = ct.DEGs.forw$Organ.A,
  Age = ct.DEGs.forw$Age,
  N = -ct.DEGs.forw$N
)


ct.DEGs.identi <- data.table( Organ.A=rep(organs,each=length(ages)),
                              Organ.B=rep(organs,each=length(ages)),
                              Age=rep(ages,times=length(organs)),
                              N=0 )
ct.DEGs <- rbindlist(list(ct.DEGs.forw,ct.DEGs.rev,ct.DEGs.identi))

rm(ct.DEGs.Up.forw,ct.DEGs.Down.forw,ct.DEGs.forw,ct.DEGs.rev,ct.DEGs.identi)




# plot circos -------------------------------------------------------------

ct.DEGs$Organ.A <- factor(ct.DEGs$Organ.A,levels=organs)
ct.DEGs$Organ.B <- factor(ct.DEGs$Organ.B,levels=organs)

# Initialize
circos.clear()
pdf("./charts/Fig.2a.CircosPlot.cpm.rev.pdf",width=15,height=15)
par(mar = c(5, 5, 5, 5), lwd = 0.1, cex=0.6)

circos.par(gap.after = c(2,2,2,2,2,2,2,2,2,2,12))
circos.initialize(factors = ct.DEGs$Organ.B, xlim=c(0,length(organs)))


# # Track A
circos.trackPlotRegion(ylim=c(0,1),track.height = 0.06, bg.col = colors.organ, bg.border='gray80',
                       panel.fun = function(x,y){
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.text(mean(xlim),ylim[2]+1.2,
                                     cex=6, 
                                     labels=get.cell.meta.data("sector.index"))})


#
# Track B - E
sapply(rev(as.character(ages)), function(currAge){
  scale.max <- max(ct.DEGs[Age==currAge]$N)*0.7
  scale.min <- min(ct.DEGs[Age==currAge]$N)*0.7
  circos.trackPlotRegion(ylim = c(scale.max, scale.min), track.height=0.15, bg.border='gray80',
                         factors = ct.DEGs[Age==currAge]$Organ.B,
                         x = as.numeric(ct.DEGs[Age==currAge]$Organ.A),
                         y = ct.DEGs[Age==currAge]$N,
                         panel.fun = function(x, y) {
                           for( i in 1:length(x)){
                             if( y[i]>0 ){
                               circos.rect(xleft=x[i]-1,xright=x[i],
                                           ybottom=0,ytop=y[i],
                                           col=colors.regulateType['Down'],border='white')
                             }
                             if( y[i]<0 ){
                               circos.rect(xleft=x[i]-1,xright=x[i],
                                           ybottom=y[i],ytop=0,
                                           col=colors.regulateType['Up'],border='white')
                             }
                             if( y[i]==0 ){
                               circos.rect(xleft=x[i]-1,xright=x[i],
                                           ybottom=-10,ytop=10,
                                           col='gray20',border='white')
                             }
                           }
                         }
  )
  return(0)
})

# Track F
circos.trackPlotRegion(ylim = c(0, 1), factors = ct.DEGs$Organ.B, track.height=0.06, bg.border='gray80',
                       x = as.numeric(ct.DEGs$Organ.A),
                       panel.fun = function(x, y) {
                         for( i in 1:length(x)){
                           circos.rect(xleft=x[i]-1,xright=x[i],ybottom=0,ytop=1,col=colors.organ[x[i]])
                         }
                       })



dev.off()



# numbers in text ---------------------------------------------------------

mean(ct.DEGs[Organ.B=='Brn'][N<0]$N)
mean(ct.DEGs[Organ.B=='Brn'][N>0]$N)

mean(ct.DEGs[Organ.B=='Utr'][N<0]$N)
mean(ct.DEGs[Organ.B=='Utr'][N>0]$N)

mean(ct.DEGs[Organ.B=='Msc'][N<0]$N)
mean(ct.DEGs[Organ.B=='Msc'][N>0]$N)

mean(ct.DEGs[Organ.B=='Tst'][N<0][Age%in%c('2', '6','21')]$N)
mean(ct.DEGs[Organ.B=='Tst'][N>0][Age%in%c('2', '6','21')]$N)
mean(ct.DEGs[Organ.B=='Tst'][N<0][Age%in%c('104')]$N)
mean(ct.DEGs[Organ.B=='Tst'][N>0][Age%in%c('104')]$N)

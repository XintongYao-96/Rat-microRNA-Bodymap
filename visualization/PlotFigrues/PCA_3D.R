
setwd("~/Downloads/Rat Bodymap/")
PCAresult <- readRDS("./data/meta_PCS.rds")
# Dependent libraries
library(data.table)
library(RColorBrewer)
library(ggplot2)
dt.meta <- readRDS('./data/Rat_SciData/data/metadata.rds')
# Set colors
organs <- unique(dt.meta$Organ)
colors.organ <- colorRampPalette(brewer.pal(9,'Set1'))(11)
names(colors.organ) <- organs

age <- unique(dt.meta$Age)
colors.age <- c(15,16,17,18)
names(colors.age) <- age

colors.organ.pastel <- colorRampPalette(brewer.pal(9,'Set1'))(11)
names(colors.organ.pastel) <- organs
colors.organ.pastel_table <- as.data.frame(colors.organ.pastel)
colors.age.pastel_table <- data.frame(Age=unique(PCAresult$Age),colors=c(16,17,15,3))
colors.age.pastel_table <- data.frame(shape=colors.age.pastel_table[,-1])
rownames(colors.age.pastel_table) <- unique(PCAresult$Age)
colnames(colors.age.pastel_table) <- "shape"
xmax <- ceiling(max(PCAresult$PC1)/10)*10
xmin <- floor(min(PCAresult$PC1)/10)*10
ymax <- ceiling(max(PCAresult$PC2)/10)*10
ymin <- floor(min(PCAresult$PC2)/10)*10
zmax <- ceiling(max(PCAresult$PC3)/10)*10
zmin <- floor(min(PCAresult$PC3)/10)*10

PCAresult$Group2 <- paste(PCAresult$Organ,PCAresult$Age,sep = "_")
list <- unique(PCAresult$Group2)
PCAresult_type <- PCAresult[which(PCAresult$Group2==list[1]),]
x <- PCAresult_type$PC1
y <- PCAresult_type$PC2
z <- PCAresult_type$PC3
Organ_1 <- as.character(unique(sapply(strsplit(as.character(PCAresult_type$Group2),"_"),function(x){paste(x[1])})))
Age_1 <- as.character(unique(sapply(strsplit(as.character(PCAresult_type$Group2),"_"),function(x){paste(x[2])})))
color_1 <- colors.organ.pastel_table[Organ_1,]
shape_1 <- colors.age.pastel_table[Age_1,]
  
par(mfrow = c(1, 1))
scatter3D(x, y, z,xlab = "PC1", ylab ="PC2", zlab = "PC3",  bty = "b2",
          xlim = c(xmin,xmax), ylim = c(ymin,ymax), zlim = c(zmin,zmax),ticktype = "detailed",
         col = color_1, pch = shape_1,groups = PCAresult_type$Organ,sphere.size=10)

for (i in 2:length(list)) {
  PCAresult_type <- PCAresult[which(PCAresult$Group2==list[i]),]
  x <- PCAresult_type$PC1
  y <- PCAresult_type$PC2
  z <- PCAresult_type$PC3
  Organ_1 <- as.character(unique(sapply(strsplit(as.character(PCAresult_type$Group2),"_"),function(x){paste(x[1])})))
  Age_1 <- as.character(unique(sapply(strsplit(as.character(PCAresult_type$Group2),"_"),function(x){paste(x[2])})))
  color_1 <- colors.organ.pastel_table[Organ_1,]
  shape_1 <- colors.age.pastel_table[Age_1,]
  
  par(mfrow = c(1, 1))
  # scatter3D(x, y, z,xlab = "PC1", ylab ="PC2", zlab = "PC3",  #bty = "g",
  #           xlim = c(xmin,xmax), ylim = c(ymin,ymax), zlim = c(zmin,zmax),
  #           col = color_1, pch = shape_1,add = TRUE)
  scatter3D(x, y, z,xlab = "PC1", ylab ="PC2", zlab = "PC3",  bty = "b2",
            xlim = c(xmin,xmax), ylim = c(ymin,ymax), zlim = c(zmin,zmax),ticktype = "detailed",
            col = color_1, pch = shape_1,groups = PCAresult_type$Organ,add = TRUE)
  print(i)
}
# 
# 
# 
# PCAresult <- as.data.frame(PCAresult)
# PCAresult$Organ <- as.character(PCAresult$Organ,levels=c(as.character(rownames(as.data.frame(colors.organ.pastel)))))
# PCAresult$Age <- as.character(PCAresult$Age,levels=c(as.character(rownames(as.data.frame(colors.age)))))
# 
# ggplot(PCAresult, aes(x=PC1, y=PC2))+
#   geom_point(aes(shape=Age,color=Organ),size=2.5)+
#   scale_shape_manual(values=c(c((as.data.frame(colors.age))[,1]))) + 
#   scale_color_manual(values=c(c((as.data.frame(colors.organ.pastel))[,1])))

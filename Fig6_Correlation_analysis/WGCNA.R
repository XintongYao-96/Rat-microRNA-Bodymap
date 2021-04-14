

library(WGCNA)  
library(magrittr)
library(data.table)
options(stringsAsFactors = FALSE);


# prepare exprmat and metadata --------------------------------------------

dt.meta <- readRDS('./RDS/metadata.rds')
dt.miRNA.exprmat <- readRDS('./RDS/exprMat/exprmat_Z_miRNA_log2CPM.rds')
dt.mRNA.exprmay <- readRDS('./RDS/exprMat/exprmat_Z_mRNA_log2FPKM.rds')

# 318 samples metadata and miRNA-mRNA exprmat
dt.metadata <- dt.meta[ID.sample%in%colnames(dt.miRNA.exprmat)]%>%data.frame(.)
rownames(dt.metadata) <- dt.metadata$ID.sample
dt.combine.expmrmat <- rbind(dt.miRNA.exprmat, dt.mRNA.exprmay)

datExpr <- t(dt.combine.expmrmat)


#  select a proper power --------------------------------------------------


powers = c(c(1:10), seq(from = 12, to=50, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;


# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



# one step network construction and module detect --------

# script from tutorial
################################
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "318RatTOM",
                       verbose = 5
                       )

table(net$colors)

# script from Yu
#We choose the power 6, suibian
##################################
# 
# softPower = 6 
# adjacency = adjacency(datExpr, power = softPower)
# # Turn adjacency into topological overlap
# TOM = TOMsimilarity(adjacency)
# dissTOM = 1-TOM
# # Call the hierarchical clustering function
# geneTree = flashClust(as.dist(dissTOM), method = "average");
# # Plot the resulting clustering tree (dendrogram)
# sizeGrWindow(12,9)
# plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#      labels = FALSE, hang = 0.04);
# 
# 



# visulization ------------------------------------------------------------


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file ='./RData/miRNA-mRNA-perTissueZ-networkConstruction-auto.RData')


library(data.table)
library(edgeR)
library(limma)

Func.CallDEGs <- function(exprMat, groupInfo){
  design <- model.matrix(~Group, data=groupInfo)
  dge <- DGEList(exprMat)
  
  dge <- calcNormFactors(dge,method='none')
  v <- voom(dge, design=design)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  
  result <- suppressMessages(topTable(fit, n=Inf, adjust.method = "none",sort.by = 'logFC'))

  return(result)
}



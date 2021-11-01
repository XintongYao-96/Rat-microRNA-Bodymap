library(data.table)
library(lme4)


Function.calc.PVCA <- function(exprMat, expDesign, pct_threshold = .5876, min.PC.use=3){
  # This function is for PVCA analysis of a given expression matrix with experiental design
  
  # Parameters:
  #   exprMat: a matrix of the expression profile which a row for a gene and a column for a sample 
  #   expDesign: a data frame that describe the treatment information of each sample
  #   pct_threshold:  the amount of variability desired to be explained by the principal components.  
  #                   default .5876, which is set to match the results in book chapter and SAS code. 
  #   min.PC.use: the minimum number of PCs used for variance analysis 

  # Code:
 
  # Principle Component Analysis (PCA) --------------------------------------
  
  pca_prcomp <- prcomp(t(exprMat),scale=T) # Using  singular value decomposition instead of eigen on the covariance matrix
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$ID.sample=rownames(pcs)
  
  dt.varProp.pcs <- data.table(PCX=1:nrow(pcs),
                               Percent=summary(pca_prcomp)$importance[2,],
                               AccumPercent=summary(pca_prcomp)$importance[3,])

  # To determine Components used for VCA ------------------------------------
  N.pc.use <- max(min.PC.use , min(dt.varProp.pcs[AccumPercent>pct_threshold]$PCX))


  # Variance component analysis (VCA) ---------------------------------------

  # Generate formula for LMM
  # Two-way interaction considered
  vars.Treatment <- setdiff(colnames(expDesign),'ID.sample')
  factor.oneWay <- paste("(1|",vars.Treatment,")",sep='')
  factor.twoWay <- paste("(1|",apply(combn(vars.Treatment,m=2),2,function(x)paste(x,collapse=':')),")",sep='')

  formula.right <- paste(c(factor.oneWay,factor.twoWay),collapse=" + ")
  formula <- paste("PCX ~ ",formula.right,sep='')
  
  # Calc variance 
  dt.randomEffects <- rbindlist(lapply(1:N.pc.use, function(x){
    dt.pcX.Annot <- merge(expDesign,
                   pcs[,c('ID.sample',sprintf("PC%d",x))],
                   by='ID.sample')
    setnames(dt.pcX.Annot,sprintf("PC%d",x),'PCX')
    
    Rm1ML <- lmer( formula ,
                   data=dt.pcX.Annot, 
                   REML = TRUE, verbose = FALSE, na.action = na.omit)
    randomEffects<-c(sqrt(abs(unlist(VarCorr(Rm1ML)))),attr(VarCorr(Rm1ML), "sc"))
    names(randomEffects)[length(randomEffects)]<-'Residual'
    dt.randomEffects.single <- data.table(PCX=x, 
                                          Stage=names(randomEffects),
                                          Effect=randomEffects)
    sum.effect <- sum(dt.randomEffects.single$Effect)
    dt.randomEffects.single[,Std.Effect.weighted:=Effect/sum.effect*dt.varProp.pcs[x]$Percent]
    return(dt.randomEffects.single)
  }))

  collapse.randomEffects <- dt.randomEffects[,sum(Std.Effect.weighted),by=Stage]
  setnames(collapse.randomEffects ,2,'Effect')
  
  sum.randomEffects <- sum(collapse.randomEffects$Effect)
  collapse.randomEffects[,Perc.effect:=Effect/sum.randomEffects*100]

  return(collapse.randomEffects)
}

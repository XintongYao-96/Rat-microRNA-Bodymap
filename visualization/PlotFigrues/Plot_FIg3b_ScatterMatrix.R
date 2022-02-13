rm(list = ls())
library(data.table)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(GGally)

dt.meta <- readRDS('./RDS/metadata.rds')
exprmat <- readRDS('./RDS/exprMat/exprMat_logCPM_r604c318.rds')


lst.sample.forPlot <- dt.meta[Organ=='Adr' | Organ=='Brn' ][Age==6][Sex=='F']$Colnames.miRNA
exprmat.forPlot <- exprmat[, lst.sample.forPlot]



# plot scatter matrix -----------------------------------------------------


max = as.numeric(exprmat.forPlot) %>% na.omit() %>% max()
min = as.numeric(exprmat.forPlot) %>% na.omit() %>% min()

fix_scale <- c(min, max)

panel.points <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_abline(intercept = 0 , slope = 1 , color='red4', linetype='dashed', size=0.5)+
    geom_point(alpha=0.6,size=0.3,color='DodgerBlue4')+
    scale_x_continuous(limits=fix_scale)+
    scale_y_continuous(limits=fix_scale)
}


panel.cor <- function(data,mapping){
  eval_data_col <- function(data, aes_col) {
    rlang::eval_tidy(aes_col, data)
  }
  
  
  dt <- data.table(
    Batch1 = eval_data_col(data,mapping$x),
    Batch2 = eval_data_col(data,mapping$y)
  )
  dt = na.omit(dt)
  corr = cor(dt$Batch1, dt$Batch2)
  
  #xData <- eval_data_col(data,mapping$x)
  #yData <- eval_data_col(data,mapping$y)
  
  # corr = cor(xData,yData)
  corr.text <-  round(corr,2)
  corr.color <- brewer.pal(9,'Blues')[8]
  #corr.color <- '#0072B2'
  corr.size <- 6
  
  ggplot(data=data,mapping=mapping)+
    geom_text(x=0.5,y=0.5,label=corr.text,size=corr.size,color=corr.color)
}


ggpairs(as.data.frame(exprmat.forPlot),
        lower = list(continuous = panel.points),
        diag = list(continuous='barDiag'),
        upper = list(continuous = panel.cor))+
  theme_classic()+
  theme(panel.border = element_rect(linetype = "dashed", colour = "black", fill = NA))
ggsave('./charts/Fig3b_Scattermatrix_correlation.pdf', width = 8, height = 8)

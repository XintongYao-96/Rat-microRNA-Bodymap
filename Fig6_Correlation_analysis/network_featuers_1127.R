library(data.table)
library(magrittr)

dt.network <- fread('./supplementary/Fig6B_cluster1_miRNA_gene_pathway.csv')
dt.network.cluster1 <- dt.network[selected=='TRUE']
dt.network.cluster2 <- dt.network[selected=='FALSE']

lst.network1.miRNA <- dt.network.cluster1[type=='miRNA']$name
lst.network1.gene <- dt.network.cluster1[type=='gene']$name
lst.network1.pathway <- dt.network.cluster1[type=='pathway']$name

lst.network2.miRNA <- dt.network.cluster2[type=='miRNA']$name
lst.network2.gene <- dt.network.cluster2[type=='gene']$name
lst.network2.pathway <- dt.network.cluster2[type=='pathway']$name

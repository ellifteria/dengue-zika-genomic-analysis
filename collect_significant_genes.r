significant_gene_data <- read.csv("ZIKV_infected_results_significant.csv")

write.table(subset(significant_gene_data, log2FoldChange > 0)[1], "ZIKV_infected_upreg_gene_names.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(subset(significant_gene_data, log2FoldChange < 0)[1], "ZIKV_infected_downreg_gene_names.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

significant_gene_data <- read.csv("DENV_infected_results_significant.csv")

write.table(subset(significant_gene_data, log2FoldChange > 0)[1], "DENV_infected_upreg_gene_names.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(subset(significant_gene_data, log2FoldChange < 0)[1], "DENV_infected_downreg_gene_names.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

significant_gene_data <- read.csv("./deseq/zikv-deseq/ZIKV_infected_results_significant.csv")

write.table(subset(significant_gene_data, log2FoldChange > 0)[1], "./gene-ontology/ZIKV_infected_upreg_gene_names.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(subset(significant_gene_data, log2FoldChange < 0)[1], "./gene-ontology/ZIKV_infected_downreg_gene_names.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

significant_gene_data <- read.csv("./deseq/denv-deseq/DENV_infected_results_significant.csv")

write.table(subset(significant_gene_data, log2FoldChange > 0)[1], "./gene-ontology/DENV_infected_upreg_gene_names.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(subset(significant_gene_data, log2FoldChange < 0)[1], "./gene-ontology/DENV_infected_downreg_gene_names.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

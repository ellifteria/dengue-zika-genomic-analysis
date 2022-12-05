significant_gene_data <- read.csv("DENV_infected_results_significant_high_p.csv")

write.table(subset(significant_gene_data, log2FoldChange < 0)[1], "./gene-ontology/DENV_infected_downreg_gene_names_high_p.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
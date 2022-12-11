library("pheatmap")
range <- 1
dat <- read.table("./gene-tpms/final_gene_tpm_mod.csv", header = TRUE, sep= ",")
donor_cols <- data.frame(donor = rep(c("B", "A"), 3), infection = rep(c("DENV", "ZIKV", "uninfected"), c(2,2,2)))
rownames(donor_cols) <- colnames(dat)
pheatmap(dat, breaks=seq(-range, range, length.out = 100), fontsize=5, fontsize_row = 2.45, annotation_col = donor_cols, show_rownames = TRUE, show_colnames = FALSE, treeheight_row = 25, treeheight_col = 20)

if(file.exists("./Rplots.pdf")){
	file.rename(from="./Rplots.pdf", to="./data-visualization/heatmaps/Rplots.pdf")
}
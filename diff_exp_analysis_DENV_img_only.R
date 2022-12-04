library("DESeq2") # load in DESeq 2 library

gene.counts.file = "./gene_counts.csv"
gene.annotations.file = "./gene_annotations copy.csv"

cts <- as.matrix(read.csv(gene.counts.file,
			   row.names = "gene_id"
			  )) # restructure gene count data into an R matrix


coldata <- read.csv(gene.annotations.file, row.names=1) # read in column data from the gene annotation data
coldata <- coldata[, c("donor", "infection_status")] # restructure coldata into donor and infection_status
coldata$donor <- factor(coldata$donor) # factor column data by donor
coldata$infection_status <- factor(coldata$infection_status) # factor column data by infection_status
rownames(coldata) <- sub("fb", "", rownames(coldata)) # assign row names for column data

print(rownames(coldata))

cts <- cts[, rownames(coldata)] # add to gene count matrix

cts[1:5, ] # display first 5 rows of gene count data

coldata # display all column data

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~infection_status) # prepare gene count matrix for DESeq analysis

dds$infection_status <- relevel(dds$infection_status, ref = "uninfected")

dds <- DESeq(dds)  # run DESeq analysis on data matrix
res <- results(dds) # generate a results table for the DESeq analysis performed on the gene count data
res  # display the results table for the DESeq analysis

res <- results(dds, contrast=c("infection_status", "DENV_infected", "uninfected")) # generate a results table for the DESeq analysis with the contrast specified
res # display the results table

resultsNames(dds) # print the coefficients of the results table

resLFC <- lfcShrink(dds, coef="infection_status_DENV_infected_vs_uninfected", type="apeglm") # shrink the effect size for the coefficient `donor_B_vs_A` to facilitate easier visualization and ranking
resLFC # display the results table after effect size shrinkage

resOrdered <- res[order(res$pvalue),] # order the results by smallest p-value
summary(res) # display a summary of the results table

sum(res$padj < 0.2, na.rm=TRUE) # calculate the number of adjusted p-values less than 0.1

res05 <- results(dds, alpha=0.05) # generate a new results table with alpha set to 0.05 instead of the default value of 0.1
summary(res05) # display a summary of the new results table

sum(res05$padj < 0.05, na.rm=TRUE) # calculate the number of adjusted p-values less than 0.05

library("IHW") # load the IHW (independent hypothesis weighting) library

resIHW <- results(dds, filterFun=ihw) # generate a results table filtered by the independent hypothesis weighting function
summary(resIHW) # display a summary of the results table

sum(resIHW$padj < 0.1, na.rm=TRUE) # calculate the number of adjusted p-values less than 0.1

metadata(resIHW)$ihwResult # display the metadata for independent hypothesis weighting results of the results table 

plotMA(resLFC, ylim=c(-2, 2), main="Log2 fold change of normalized gene counts:\nDENV infected vs. uninfected") # plot the log2 fold change for the mean of the normalized counts of all the samples in the gene count DESeq analysis results table post-shrinkage

resultsNames(dds) # display the coefficients of the gene count DESeq analysis

resOrdered # display the ordered results

resOrdered[c("ENSG00000172216", "ENSG00000175197"),]